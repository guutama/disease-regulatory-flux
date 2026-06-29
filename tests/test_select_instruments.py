"""Unit tests for bin/select_instruments.py (instrument selection / BioFindr inputs).

A small simulated trio is generated in tmp_path. We test:
  * lead-SNP selection and the deterministic tie-break (min p -> largest |beta| -> variant),
  * the dE instrument list (one row per regulator, variant before gene),
  * dG: only the selected instruments, transposed to samples x SNPs, deduplicated,
  * dX: all genes, transposed to samples x genes,
  * that dG sample rows are reordered to match the expression sample order,
  * value pass-through (no transformation),
  * the QC counts,
  * the error paths (empty eQTL, sample mismatch, instrument missing from genotype).
"""
import gzip
import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "select_instruments", Path(__file__).resolve().parents[1] / "bin" / "select_instruments.py")
si = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(si)


# ----------------------------------------------------------------------------- fixtures
def _write(path: Path, rows: list[list[str]]) -> None:
    text = "\n".join("\t".join(map(str, r)) for r in rows) + "\n"
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def make_trio(d: Path):
    """Expression (4 genes x 3 samples), genotype (3 variants x 3 samples, in a DIFFERENT
    sample order to test reordering), eQTL designed so that: geneA->rs1 (clear min p),
    geneB->rs2, geneC ties on p (rs3 vs rs1) and |beta| picks rs1, geneD ties on p and
    |beta| so the first variant id (rs2) wins. rs1 is shared by geneA and geneC, so dG must
    deduplicate to {rs1, rs2}; rs3 is never a lead and is dropped from dG."""
    expr = [["gene_id", "S1", "S2", "S3"],
            ["geneA", 0.1, 0.2, 0.3],
            ["geneB", 1.1, 1.2, 1.3],
            ["geneC", 2.1, 2.2, 2.3],
            ["geneD", 3.1, 3.2, 3.3]]
    geno = [["variant_id", "S2", "S1", "S3"],     # sample order differs from expression
            ["rs1", 1, 0, 2],                     # -> expression order S1,S2,S3 = 0,1,2
            ["rs2", 1, 1, 0],                     # -> 1,1,0
            ["rs3", 9, 9, 9]]                     # not a lead -> dropped
    eqtl = [["gene_id", "variant_id", "beta", "se", "pvalue"],
            ["geneA", "rs1", 0.5, 0.1, 1e-8],
            ["geneA", "rs2", -0.2, 0.1, 1e-3],
            ["geneB", "rs2", 0.3, 0.1, 1e-5],
            ["geneC", "rs3", 0.4, 0.1, 1e-6],
            ["geneC", "rs1", -0.6, 0.1, 1e-6],    # tie on p with rs3; |beta| 0.6 > 0.4 -> rs1
            ["geneD", "rs2", 0.2, 0.1, 1e-4],
            ["geneD", "rs3", 0.2, 0.1, 1e-4]]     # tie on p and |beta|; first variant rs2 wins
    _write(d / "expr.tsv", expr)
    _write(d / "geno.tsv", geno)
    _write(d / "eqtl.tsv", eqtl)
    return d / "expr.tsv", d / "geno.tsv", d / "eqtl.tsv"


def run(d: Path, expr, geno, eqtl, extra=None):
    outs = {k: d / f"{k}" for k in ("dx.csv", "dg.csv", "de.csv", "qc.tsv")}
    argv = ["--expression", str(expr), "--genotype", str(geno), "--eqtl", str(eqtl),
            "--out-dx", str(outs["dx.csv"]), "--out-dg", str(outs["dg.csv"]),
            "--out-de", str(outs["de.csv"]), "--qc", str(outs["qc.tsv"])]
    if extra:
        argv += extra
    si.main(argv)
    return outs


def read_csv(path: Path):
    return [ln.split(",") for ln in path.read_text().splitlines()]


def read_qc(path: Path):
    return {k: v for k, v in (ln.split("\t") for ln in path.read_text().splitlines()[1:])}


# -------------------------------------------------------------------------------- tests
def test_select_leads_tiebreak():
    # rows as (gene, variant, beta, se, pvalue); indices 0,1,2,_,4
    rows = [["geneC", "rs3", "0.4", "0.1", "1e-6"],
            ["geneC", "rs1", "-0.6", "0.1", "1e-6"],   # same p, larger |beta| -> rs1
            ["geneD", "rs2", "0.2", "0.1", "1e-4"],
            ["geneD", "rs3", "0.2", "0.1", "1e-4"]]     # same p and |beta| -> first variant rs2
    leads, n_ties = si.select_leads(rows, gene_i=0, var_i=1, beta_i=2, p_i=4)
    assert leads == {"geneC": "rs1", "geneD": "rs2"}
    assert n_ties == 2


def test_instrument_list_de(tmp_path):
    expr, geno, eqtl = make_trio(tmp_path)
    outs = run(tmp_path, expr, geno, eqtl)
    de = read_csv(outs["de.csv"])
    assert de[0] == ["variant_id", "gene_id"]          # variant first, gene second
    pairs = {(r[0], r[1]) for r in de[1:]}
    assert pairs == {("rs1", "geneA"), ("rs2", "geneB"), ("rs1", "geneC"), ("rs2", "geneD")}


def test_dg_instruments_only_transposed_and_reordered(tmp_path):
    expr, geno, eqtl = make_trio(tmp_path)
    outs = run(tmp_path, expr, geno, eqtl)
    dg = read_csv(outs["dg.csv"])
    # only the selected instruments, deduplicated; rs3 dropped
    assert set(dg[0]) == {"rs1", "rs2"}
    # samples are rows, in the expression order S1,S2,S3 (genotype was stored S2,S1,S3)
    col = {name: i for i, name in enumerate(dg[0])}
    rows = dg[1:]
    assert [r[col["rs1"]] for r in rows] == ["0", "1", "2"]
    assert [r[col["rs2"]] for r in rows] == ["1", "1", "0"]


def test_dx_all_genes_transposed(tmp_path):
    expr, geno, eqtl = make_trio(tmp_path)
    outs = run(tmp_path, expr, geno, eqtl)
    dx = read_csv(outs["dx.csv"])
    assert dx[0] == ["geneA", "geneB", "geneC", "geneD"]   # all genes are columns
    assert len(dx) == 1 + 3                                # header + 3 samples
    col = {name: i for i, name in enumerate(dx[0])}
    # first sample row S1 carries each gene's S1 value, unchanged
    assert dx[1][col["geneA"]] == "0.1"
    assert dx[1][col["geneD"]] == "3.1"
    # last sample row S3
    assert dx[3][col["geneC"]] == "2.3"


def test_qc_counts(tmp_path):
    expr, geno, eqtl = make_trio(tmp_path)
    outs = run(tmp_path, expr, geno, eqtl)
    qc = read_qc(outs["qc.tsv"])
    assert qc["samples"] == "3"
    assert qc["genes_dX"] == "4"
    assert qc["regulators_dE"] == "4"
    assert qc["instruments_dG"] == "2"      # rs1, rs2 (deduplicated)
    assert qc["pvalue_ties_broken"] == "2"


def test_gzip_roundtrip(tmp_path):
    expr, geno, eqtl = make_trio(tmp_path)
    gz = tmp_path / "gz"; gz.mkdir()
    for name in ("expr", "geno", "eqtl"):
        (gz / f"{name}.tsv.gz").write_bytes(gzip.compress((tmp_path / f"{name}.tsv").read_bytes()))
    dx, dg, de, qc = (gz / "dx.csv.gz", gz / "dg.csv.gz", gz / "de.csv.gz", gz / "qc.tsv")
    si.main(["--expression", str(gz / "expr.tsv.gz"), "--genotype", str(gz / "geno.tsv.gz"),
             "--eqtl", str(gz / "eqtl.tsv.gz"), "--out-dx", str(dx), "--out-dg", str(dg),
             "--out-de", str(de), "--qc", str(qc)])
    header = gzip.decompress(dx.read_bytes()).decode().splitlines()[0]
    assert header.split(",") == ["geneA", "geneB", "geneC", "geneD"]


def test_empty_eqtl_errors(tmp_path):
    expr, geno, _ = make_trio(tmp_path)
    _write(tmp_path / "empty.tsv", [["gene_id", "variant_id", "beta", "se", "pvalue"]])
    with pytest.raises(SystemExit, match="no rows"):
        run(tmp_path, expr, geno, tmp_path / "empty.tsv")


def test_sample_mismatch_errors(tmp_path):
    expr = [["gene_id", "A1", "A2"], ["geneA", 1, 2]]
    geno = [["variant_id", "B1", "B2"], ["rs1", 0, 1]]
    eqtl = [["gene_id", "variant_id", "beta", "se", "pvalue"], ["geneA", "rs1", 0.1, 0.1, 1e-3]]
    _write(tmp_path / "e.tsv", expr); _write(tmp_path / "g.tsv", geno); _write(tmp_path / "q.tsv", eqtl)
    with pytest.raises(SystemExit, match="same samples"):
        run(tmp_path, tmp_path / "e.tsv", tmp_path / "g.tsv", tmp_path / "q.tsv")


def test_instrument_missing_from_genotype_errors(tmp_path):
    expr = [["gene_id", "S1", "S2"], ["geneA", 1, 2]]
    geno = [["variant_id", "S1", "S2"], ["rsOTHER", 0, 1]]
    eqtl = [["gene_id", "variant_id", "beta", "se", "pvalue"], ["geneA", "rs1", 0.1, 0.1, 1e-3]]
    _write(tmp_path / "e.tsv", expr); _write(tmp_path / "g.tsv", geno); _write(tmp_path / "q.tsv", eqtl)
    with pytest.raises(SystemExit, match="absent from the genotype"):
        run(tmp_path, tmp_path / "e.tsv", tmp_path / "g.tsv", tmp_path / "q.tsv")
