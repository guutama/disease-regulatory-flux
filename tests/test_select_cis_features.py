"""Unit tests for bin/select_cis_features.py (cis-eQTL feature selection).

A small eQTL table is generated in tmp_path. We test:
  * the per-gene cis feature set and that duplicate gene-variant pairs are dropped,
  * the optional p-value filter when a p-value column is present,
  * that the threshold is ignored (every pair trusted) when there is no p-value column,
  * the per-gene summary and QC counts,
  * gzip round-trip and the empty-input error.
"""
import gzip
import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "select_cis_features", Path(__file__).resolve().parents[1] / "bin" / "select_cis_features.py")
cf = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(cf)


def _write(path: Path, rows: list[list[str]]) -> None:
    text = "\n".join("\t".join(map(str, r)) for r in rows) + "\n"
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def run(d: Path, eqtl: Path, extra=None):
    argv = ["--eqtl", str(eqtl), "--tissue", "T", "--outdir", str(d)]
    if extra:
        argv += extra
    cf.main(argv)
    return {
        "features": d / "T.cis_features.tsv.gz",
        "summary": d / "T.cis_features_summary.tsv",
        "qc": d / "T.cis_features_qc.tsv",
    }


def read_tsv_gz(path: Path):
    return [ln.split("\t") for ln in gzip.decompress(path.read_bytes()).decode().splitlines()]


def read_tsv(path: Path):
    return [ln.split("\t") for ln in path.read_text().splitlines()]


def read_qc(path: Path):
    return {k: v for k, v in (ln.split("\t") for ln in path.read_text().splitlines()[1:])}


# eQTL with a p-value column: geneA has a duplicate pair and one weak (high-p) pair.
EQTL_P = [["gene_id", "variant_id", "beta", "se", "pvalue"],
          ["geneA", "rs1", 0.5, 0.1, 1e-8],
          ["geneA", "rs2", -0.2, 0.1, 1e-2],
          ["geneA", "rs1", 0.5, 0.1, 1e-8],     # duplicate -> dropped
          ["geneB", "rs3", 0.3, 0.1, 1e-6],
          ["geneB", "rs4", 0.4, 0.1, 0.5]]       # weak


def test_feature_set_and_dedup(tmp_path):
    eqtl = tmp_path / "eqtl.tsv"; _write(eqtl, EQTL_P)
    outs = run(tmp_path, eqtl)
    pairs = [(r[0], r[1]) for r in read_tsv_gz(outs["features"])[1:]]
    assert pairs == [("geneA", "rs1"), ("geneA", "rs2"), ("geneB", "rs3"), ("geneB", "rs4")]
    qc = read_qc(outs["qc"])
    assert qc["pairs_in"] == "5"
    assert qc["pairs_kept"] == "4"
    assert qc["dropped_duplicate"] == "1"
    assert qc["pvalue_filter_applied"] == "no"


def test_pvalue_filter(tmp_path):
    eqtl = tmp_path / "eqtl.tsv"; _write(eqtl, EQTL_P)
    outs = run(tmp_path, eqtl, ["--max-pvalue", "1e-3"])
    pairs = {(r[0], r[1]) for r in read_tsv_gz(outs["features"])[1:]}
    # rs2 (1e-2) and rs4 (0.5) drop; rs1 (1e-8) and rs3 (1e-6) stay
    assert pairs == {("geneA", "rs1"), ("geneB", "rs3")}
    qc = read_qc(outs["qc"])
    assert qc["dropped_pvalue"] == "2"
    assert qc["pvalue_filter_applied"] == "yes"


def test_pvalue_threshold_ignored_without_pvalue_column(tmp_path):
    eqtl = tmp_path / "eqtl.tsv"
    _write(eqtl, [["gene_id", "variant_id"], ["geneA", "rs1"], ["geneB", "rs3"]])
    outs = run(tmp_path, eqtl, ["--max-pvalue", "1e-3"])    # no p-value col -> trust all
    pairs = {(r[0], r[1]) for r in read_tsv_gz(outs["features"])[1:]}
    assert pairs == {("geneA", "rs1"), ("geneB", "rs3")}
    assert read_qc(outs["qc"])["pvalue_filter_applied"] == "no"


def test_summary_counts(tmp_path):
    eqtl = tmp_path / "eqtl.tsv"; _write(eqtl, EQTL_P)
    outs = run(tmp_path, eqtl)
    summ = {r[0]: r[1] for r in read_tsv(outs["summary"])[1:]}
    assert summ == {"geneA": "2", "geneB": "2"}


def test_gzip_input(tmp_path):
    eqtl = tmp_path / "eqtl.tsv.gz"; _write(eqtl, EQTL_P)
    outs = run(tmp_path, eqtl)
    assert read_tsv_gz(outs["features"])[0] == ["gene_id", "variant_id"]


def test_empty_input_errors(tmp_path):
    eqtl = tmp_path / "eqtl.tsv"
    _write(eqtl, [["gene_id", "variant_id", "pvalue"]])
    with pytest.raises(SystemExit, match="no cis feature set|No cis-eQTL"):
        run(tmp_path, eqtl)
