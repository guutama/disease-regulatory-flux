"""Unit tests for bin/harmonise_inputs.py (Step 1: input harmonisation).

Simulated trio is generated in tmp_path — no data files are committed. We test:
  * the variant-key resolution rule (prefer rsID, else chr:pos:ref:alt, else error),
  * the three-way intersection (samples / genes / variants),
  * that a gene losing all its instruments is dropped,
  * sample-order consistency between aligned expression and genotype,
  * value pass-through (no transformation),
  * the QC counts,
  * the empty-overlap error paths.
"""
import gzip
import importlib.util
import sys
from pathlib import Path

import pytest

# import the script as a module
_SPEC = importlib.util.spec_from_file_location(
    "harmonise_inputs", Path(__file__).resolve().parents[1] / "bin" / "harmonise_inputs.py")
hi = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(hi)


# ----------------------------------------------------------------------------- fixtures
def _write(path: Path, rows: list[list[str]]) -> None:
    text = "\n".join("\t".join(map(str, r)) for r in rows) + "\n"
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def make_rsid_trio(d: Path):
    """Expression (4 genes x 4 samples), genotype (4 variants x 4 samples), eQTL keyed by
    rsID. Designed so that: S4 is genotype-only and S5 expression-only (so the shared
    samples are S1,S2,S3); geneD has no cis-eQTL (it must still be KEPT as a target); rsX
    is not in the eQTL map (so the genotype is reduced to rs1,rs2,rs3)."""
    expr = [["gene_id", "S1", "S2", "S3", "S5"],
            ["geneA", 0.1, 0.2, 0.3, 0.4],
            ["geneB", 1.1, 1.2, 1.3, 1.4],
            ["geneC", 2.1, 2.2, 2.3, 2.4],
            ["geneD", 3.1, 3.2, 3.3, 3.4]]
    geno = [["rsID", "S1", "S2", "S3", "S4"],
            ["rs1", 0, 1, 2, 0],
            ["rs2", 1, 1, 0, 2],
            ["rs3", 2, 0, 1, 1],
            ["rsX", 0, 0, 0, 0]]
    eqtl = [["gene_id", "rsID", "beta", "se", "pvalue"],
            ["geneA", "rs1", 0.5, 0.1, 1e-8],
            ["geneA", "rs2", -0.2, 0.1, 1e-3],
            ["geneB", "rs2", 0.3, 0.1, 1e-5],
            ["geneC", "rs3", 0.4, 0.1, 1e-6]]
    _write(d / "expr.tsv", expr)
    _write(d / "geno.tsv", geno)
    _write(d / "eqtl.tsv", eqtl)
    return d / "expr.tsv", d / "geno.tsv", d / "eqtl.tsv"


def run(d: Path, expr, geno, eqtl, extra=None):
    outs = {k: d / f"{k}.tsv" for k in ("oe", "og", "oq", "oqc")}
    argv = ["--expression", str(expr), "--genotype", str(geno), "--eqtl", str(eqtl),
            "--out-expression", str(outs["oe"]), "--out-genotype", str(outs["og"]),
            "--out-eqtl", str(outs["oq"]), "--qc", str(outs["oqc"])]
    if extra:
        argv += extra
    hi.main(argv)
    return outs


def read_tsv(path: Path):
    return [ln.split("\t") for ln in path.read_text().splitlines()]


def read_qc(path: Path):
    return {k: v for k, v in (ln.split("\t") for ln in path.read_text().splitlines()[1:])}


# -------------------------------------------------------------------------------- tests
def test_rsid_key_resolution():
    assert hi.resolve_key_type(["rsID", "S1"], ["gene_id", "rsID"],
                               hi.parse_args(_min_args())) == "rsid"


def test_composite_key_resolution():
    args = hi.parse_args(_min_args())
    geno = ["chromosome", "position", "ref", "alt", "S1"]
    eqtl = ["gene_id", "chromosome", "position", "ref", "alt", "beta"]
    assert hi.resolve_key_type(geno, eqtl, args) == "composite"


def test_no_shared_key_errors():
    args = hi.parse_args(_min_args())
    with pytest.raises(SystemExit):
        hi.resolve_key_type(["rsID", "S1"], ["gene_id", "chromosome", "position", "ref", "alt"], args)


def test_full_alignment(tmp_path):
    expr, geno, eqtl = make_rsid_trio(tmp_path)
    outs = run(tmp_path, expr, geno, eqtl)

    # samples: expr {S1,S2,S3,S5} ∩ geno {S1..S4} = {S1,S2,S3}, in expression order
    assert read_tsv(outs["oe"])[0] == ["gene_id", "S1", "S2", "S3"]
    assert read_tsv(outs["og"])[0] == ["variant_id", "S1", "S2", "S3"]

    # ALL genes kept, incl. geneD which has no cis-eQTL (a target, not a regulator)
    genes = [r[0] for r in read_tsv(outs["oe"])[1:]]
    assert set(genes) == {"geneA", "geneB", "geneC", "geneD"}
    # genotype reduced to the eQTL variants; rsX (not in eQTL) dropped
    variants = [r[0] for r in read_tsv(outs["og"])[1:]]
    assert set(variants) == {"rs1", "rs2", "rs3"}

    # eQTL pairs only for kept genes & variants
    pairs = {(r[0], r[1]) for r in read_tsv(outs["oq"])[1:]}
    assert pairs == {("geneA", "rs1"), ("geneA", "rs2"), ("geneB", "rs2"), ("geneC", "rs3")}


def test_value_passthrough_and_sample_reorder(tmp_path):
    expr, geno, eqtl = make_rsid_trio(tmp_path)
    outs = run(tmp_path, expr, geno, eqtl)
    # geneA dosages for S1,S2,S3 unchanged; S4 column dropped
    rs1 = next(r for r in read_tsv(outs["og"])[1:] if r[0] == "rs1")
    assert rs1[1:] == ["0", "1", "2"]
    # geneA expression for S1,S2,S3 unchanged; S5 (no genotype) dropped
    geneA = next(r for r in read_tsv(outs["oe"])[1:] if r[0] == "geneA")
    assert geneA[1:] == ["0.1", "0.2", "0.3"]


def test_qc_counts(tmp_path):
    expr, geno, eqtl = make_rsid_trio(tmp_path)
    outs = run(tmp_path, expr, geno, eqtl)
    qc = read_qc(outs["oqc"])
    assert qc["variant_key_type"] == "rsid"
    assert qc["samples_kept"] == "3"
    assert qc["genes_written"] == "4"        # all expression genes kept
    assert qc["regulators_kept"] == "3"      # geneA/B/C have a cis-eQTL; geneD does not
    assert qc["variants_kept"] == "3"
    assert qc["eqtl_pairs_kept"] == "4"
    assert qc["variants_genotype"] == "4"


def test_gzip_roundtrip(tmp_path):
    expr, geno, eqtl = make_rsid_trio(tmp_path)
    # re-write inputs gzipped, outputs gzipped
    gz = tmp_path / "gz"; gz.mkdir()
    for name in ("expr", "geno", "eqtl"):
        (gz / f"{name}.tsv.gz").write_bytes(gzip.compress((tmp_path / f"{name}.tsv").read_bytes()))
    outs = {k: gz / f"{k}.tsv.gz" for k in ("oe", "og", "oq")}
    qc = gz / "qc.tsv"
    hi.main(["--expression", str(gz / "expr.tsv.gz"), "--genotype", str(gz / "geno.tsv.gz"),
             "--eqtl", str(gz / "eqtl.tsv.gz"), "--out-expression", str(outs["oe"]),
             "--out-genotype", str(outs["og"]), "--out-eqtl", str(outs["oq"]), "--qc", str(qc)])
    rows = [ln.split("\t") for ln in gzip.decompress(outs["oe"].read_bytes()).decode().splitlines()]
    assert rows[0] == ["gene_id", "S1", "S2", "S3"]


def test_empty_sample_overlap_errors(tmp_path):
    expr = [["gene_id", "A1", "A2"], ["geneA", 1, 2]]
    geno = [["rsID", "B1", "B2"], ["rs1", 0, 1]]
    eqtl = [["gene_id", "rsID", "beta", "se", "pvalue"], ["geneA", "rs1", 0.1, 0.1, 1e-3]]
    _write(tmp_path / "e.tsv", expr); _write(tmp_path / "g.tsv", geno); _write(tmp_path / "q.tsv", eqtl)
    with pytest.raises(SystemExit, match="samples"):
        run(tmp_path, tmp_path / "e.tsv", tmp_path / "g.tsv", tmp_path / "q.tsv")


def _min_args():
    # the smallest valid argv so parse_args() succeeds (paths are unused in key-resolution tests)
    return ["--expression", "x", "--genotype", "x", "--eqtl", "x",
            "--out-expression", "x", "--out-genotype", "x", "--out-eqtl", "x", "--qc", "x"]
