"""Unit tests for bin/ld_prune.py (LD filtering of the per-gene cis-SNP set).

A small genotype with known LD structure is generated in tmp_path. We test:
  * that a SNP in perfect LD with an already-kept one is dropped,
  * that p-value ordering keeps the more significant SNP of an LD pair,
  * that a constant (zero-variance) SNP is dropped,
  * the 0/1/2 hard-call matrix: dosage rounding and restriction to the kept SNP universe,
  * the per-gene summary and QC counts,
  * the error when a cis SNP is absent from the genotype matrix.
"""
import gzip
import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "ld_prune", Path(__file__).resolve().parents[1] / "bin" / "ld_prune.py")
lp = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(lp)


def _write(path: Path, rows: list[list[str]]) -> None:
    text = "\n".join("\t".join(map(str, r)) for r in rows) + "\n"
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


S = [f"S{i}" for i in range(1, 11)]                       # 10 samples

# rs1 == rs2 (perfect LD); rs3 independent of rs1; rs4 constant; rs5 fractional dosages;
# rs6 == rs7 (perfect LD).
GENO = [["variant_id", "chromosome", "position", "ref", "alt"] + S,
        ["rs1", "1", "101", "A", "G", 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        ["rs2", "1", "102", "C", "T", 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
        ["rs3", "1", "103", "G", "A", 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
        ["rs4", "1", "104", "T", "C", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ["rs5", "2", "201", "A", "C", 0.4, 1.6, 2.0, 0.0, 1.0, 0.6, 1.4, 0.9, 0.1, 1.9],
        ["rs6", "2", "202", "G", "T", 2, 2, 2, 1, 1, 1, 0, 0, 0, 2],
        ["rs7", "2", "203", "C", "A", 2, 2, 2, 1, 1, 1, 0, 0, 0, 2]]

CIS = [["gene_id", "variant_id"],
       ["geneA", "rs1"], ["geneA", "rs2"], ["geneA", "rs3"],
       ["geneB", "rs4"], ["geneB", "rs5"],
       ["geneC", "rs6"], ["geneC", "rs7"]]

# p-values: in geneA, rs1 beats rs2 (so rs1 is kept); in geneC, rs7 beats rs6.
EQTL = [["gene_id", "variant_id", "pvalue"],
        ["geneA", "rs1", 1e-8], ["geneA", "rs2", 1e-4], ["geneA", "rs3", 1e-6],
        ["geneB", "rs4", 1e-3], ["geneB", "rs5", 1e-5],
        ["geneC", "rs6", 1e-2], ["geneC", "rs7", 1e-7]]


def run(d: Path, *, eqtl=True, r2=0.7):
    _write(d / "cis.tsv", CIS)
    _write(d / "geno.tsv", GENO)
    argv = ["--cis-features", str(d / "cis.tsv"), "--genotype", str(d / "geno.tsv"),
            "--tissue", "T", "--outdir", str(d), "--ld-r2-threshold", str(r2)]
    if eqtl:
        _write(d / "eqtl.tsv", EQTL)
        argv += ["--eqtl", str(d / "eqtl.tsv")]
    lp.main(argv)
    return {
        "pruned": d / "T.cis_snps_pruned.tsv.gz",
        "geno012": d / "T.genotype_012.tsv.gz",
        "summary": d / "T.ld_prune_summary.tsv",
        "qc": d / "T.ld_prune_qc.tsv",
    }


def read_tsv_gz(path: Path):
    return [ln.split("\t") for ln in gzip.decompress(path.read_bytes()).decode().splitlines()]


def read_tsv(path: Path):
    return [ln.split("\t") for ln in path.read_text().splitlines()]


def kept_by_gene(pruned: Path):
    out: dict[str, list[str]] = {}
    for g, v in read_tsv_gz(pruned)[1:]:
        out.setdefault(g, []).append(v)
    return out


def test_drops_perfect_ld_and_keeps_independent(tmp_path):
    kept = kept_by_gene(run(tmp_path)["pruned"])
    # rs1 kept (lowest p), rs3 kept (independent), rs2 dropped (LD with rs1)
    assert set(kept["geneA"]) == {"rs1", "rs3"}


def test_pvalue_ordering_keeps_most_significant(tmp_path):
    kept = kept_by_gene(run(tmp_path)["pruned"])
    # geneC: rs7 (p=1e-7) beats rs6 (p=1e-2); the kept one is the more significant
    assert kept["geneC"] == ["rs7"]


def test_constant_snp_dropped(tmp_path):
    outs = run(tmp_path)
    kept = kept_by_gene(outs["pruned"])
    assert kept["geneB"] == ["rs5"]            # rs4 is constant -> dropped
    summ = {r[0]: r for r in read_tsv(outs["summary"])[1:]}
    # geneB row: n_cis_in=2, n_cis_kept=1, n_dropped_ld=0, n_dropped_constant=1
    assert summ["geneB"][1:] == ["2", "1", "0", "1"]


def test_hardcall_rounding_and_universe(tmp_path):
    outs = run(tmp_path)
    geno = read_tsv_gz(outs["geno012"])
    assert geno[0] == ["variant_id", "chromosome", "position", "ref", "alt"] + S
    # allele annotation is carried through for each kept variant
    annot = {r[0]: r[1:5] for r in geno[1:]}
    assert annot["rs5"] == ["2", "201", "A", "C"]
    rows = {r[0]: r[5:] for r in geno[1:]}
    # only the kept universe is written: rs1, rs3, rs5, rs7 (rs2, rs4, rs6 dropped)
    assert set(rows) == {"rs1", "rs3", "rs5", "rs7"}
    # rs5 dosages rounded to 0/1/2 (<0.5->0, [0.5,1.5)->1, >=1.5->2)
    assert rows["rs5"] == ["0", "2", "2", "0", "1", "1", "1", "1", "0", "2"]


def test_qc_counts(tmp_path):
    outs = run(tmp_path)
    qc = {r[0]: r[1] for r in read_tsv(outs["qc"])[1:]}
    assert qc["genes"] == "3"
    assert qc["samples"] == "10"
    assert qc["cis_snps_in_universe"] == "7"
    assert qc["kept_snp_universe"] == "4"
    assert qc["ordered_by"] == "pvalue"


def test_without_eqtl_uses_input_order(tmp_path):
    outs = run(tmp_path, eqtl=False)
    qc = {r[0]: r[1] for r in read_tsv(outs["qc"])[1:]}
    assert qc["ordered_by"] == "input_order"
    # rs1 comes before rs2 in the cis-feature file, so rs1 is kept and rs2 dropped
    assert set(kept_by_gene(outs["pruned"])["geneA"]) == {"rs1", "rs3"}


def test_missing_genotype_errors(tmp_path):
    _write(tmp_path / "cis.tsv", [["gene_id", "variant_id"], ["geneA", "rsX"]])
    _write(tmp_path / "geno.tsv", [["variant_id", "chromosome", "position", "ref", "alt"] + S,
                                   ["rs1", "1", "101", "A", "G"] + [0] * 10])
    with pytest.raises(SystemExit, match="absent from the genotype"):
        lp.main(["--cis-features", str(tmp_path / "cis.tsv"),
                 "--genotype", str(tmp_path / "geno.tsv"),
                 "--tissue", "T", "--outdir", str(tmp_path)])
