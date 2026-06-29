"""Unit tests for bin/trans_features.py (trans-feature construction).

A small network and pruned cis-SNP set are generated in tmp_path. The network is
A->B, A->C, B->C, C->D (so C's parents are A and B; D's parent is C), with cis SNPs on
A, B and C and none on D. We test:
  * the per-gene trans features collected by walking up to a given hop,
  * that an ancestor reachable at two different hops is reported at each,
  * that a gene with no regulators (a root) produces no features,
  * the per-gene summary and QC counts,
  * that --max-hop limits the walk depth.
"""
import gzip
import importlib.util
from pathlib import Path

_SPEC = importlib.util.spec_from_file_location(
    "trans_features", Path(__file__).resolve().parents[1] / "bin" / "trans_features.py")
tf = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(tf)

GRN = [["Source", "Target", "Probability"],
       ["A", "B", "1.0"],
       ["A", "C", "1.0"],
       ["B", "C", "1.0"],
       ["C", "D", "1.0"]]

CIS = [["gene_id", "variant_id"],
       ["A", "a1"], ["A", "a2"],
       ["B", "b1"],
       ["C", "c1"]]            # D has no cis SNPs


def _write(path: Path, rows, delim) -> None:
    path.write_text("\n".join(delim.join(map(str, r)) for r in rows) + "\n")


def run(d: Path, max_hop=2):
    _write(d / "grn.csv", GRN, ",")
    _write(d / "cis.tsv", CIS, "\t")
    tf.main(["--grn", str(d / "grn.csv"), "--cis-pruned", str(d / "cis.tsv"),
             "--tissue", "T", "--outdir", str(d), "--max-hop", str(max_hop)])
    return {
        "features": d / "T.trans_features.tsv.gz",
        "summary": d / "T.trans_features_summary.tsv",
        "qc": d / "T.trans_features_qc.tsv",
    }


def read_features(path: Path):
    rows = [ln.split("\t") for ln in gzip.decompress(path.read_bytes()).decode().splitlines()]
    assert rows[0] == ["gene_id", "hop", "source_gene_id", "variant_id"]
    return [tuple(r) for r in rows[1:]]


def read_tsv(path: Path):
    return [ln.split("\t") for ln in path.read_text().splitlines()]


def test_walk_collects_ancestor_cis(tmp_path):
    feats = read_features(run(tmp_path)["features"])
    d_rows = {(hop, src, v) for g, hop, src, v in feats if g == "D"}
    # D: hop1 -> C (c1); hop2 -> B (b1), A (a1, a2)
    assert d_rows == {("1", "C", "c1"), ("2", "B", "b1"), ("2", "A", "a1"), ("2", "A", "a2")}


def test_ancestor_reached_at_two_hops_reported_at_each(tmp_path):
    feats = read_features(run(tmp_path)["features"])
    a_for_c = {hop for g, hop, src, v in feats if g == "C" and src == "A"}
    # A is C's direct parent (hop 1) and also reachable via B (hop 2)
    assert a_for_c == {"1", "2"}


def test_root_gene_has_no_features(tmp_path):
    feats = read_features(run(tmp_path)["features"])
    assert not any(g == "A" for g, *_ in feats)        # A has no parents
    summ_genes = {r[0] for r in read_tsv(run(tmp_path)["summary"])[1:]}
    assert "A" not in summ_genes
    assert summ_genes == {"B", "C", "D"}


def test_summary_counts(tmp_path):
    summ = {r[0]: r[1:] for r in read_tsv(run(tmp_path)["summary"])[1:]}
    # C: regulators {A, B} = 2; unique SNPs {a1, a2, b1} = 3
    assert summ["C"] == ["2", "3"]
    # D: regulators {A, B, C} = 3; unique SNPs {a1, a2, b1, c1} = 4
    assert summ["D"] == ["3", "4"]


def test_qc_counts(tmp_path):
    qc = {r[0]: r[1] for r in read_tsv(run(tmp_path)["qc"])[1:]}
    assert qc["genes_in_network"] == "4"
    assert qc["regulators_with_cis"] == "3"        # A, B, C have cis; D does not
    assert qc["genes_with_trans_features"] == "3"  # B, C, D
    assert qc["trans_feature_rows"] == "11"


def test_max_hop_limits_depth(tmp_path):
    feats = read_features(run(tmp_path, max_hop=1)["features"])
    assert {hop for _g, hop, _s, _v in feats} == {"1"}    # only direct regulators
    d_rows = {(src, v) for g, hop, src, v in feats if g == "D"}
    assert d_rows == {("C", "c1")}                         # D's hop-1 ancestor only
