"""Unit and end-to-end tests for bin/run_findr_py.py (GRN reconstruction with python findr).

The pure-Python steps -- instrument/gene pairing, q-values, sub-test combination, greedy-edges
cycle removal and the master edge writer -- are tested directly. A full end-to-end test builds a
small dataset with two instrumented regulators (A drives B, E drives D; C is null), runs the real
libfindr core, and checks the three outputs, their schemas and that findr recovers the causal
edges. The end-to-end test is skipped where the ``findr`` package / libfindr are unavailable.
"""
import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

_SPEC = importlib.util.spec_from_file_location(
    "run_findr_py", Path(__file__).resolve().parents[1] / "bin" / "run_findr_py.py")
rf = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rf)


# --------------------------------------------------------------------------- pure functions
def test_getpairs_maps_and_sorts_by_source_gene():
    genes = ["A", "B", "C"]
    snps = ["s1", "s2"]
    dE = pd.DataFrame({"variant_id": ["s2", "s1"], "gene_id": ["C", "A"]})
    pairs = rf.getpairs(genes, snps, dE)
    # (colG, colX): s2->C = (1, 2), s1->A = (0, 0); sorted by colX -> A first
    assert pairs.tolist() == [[0, 0], [1, 2]]


def test_getpairs_drops_missing_gene_and_raises_on_missing_snp():
    genes = ["A", "B"]
    snps = ["s1"]
    # gene "Z" absent from dX -> row dropped, leaving the valid s1->A pair
    ok = rf.getpairs(genes, snps, pd.DataFrame(
        {"variant_id": ["s1", "s1"], "gene_id": ["A", "Z"]}))
    assert ok.tolist() == [[0, 0]]
    # SNP "sX" absent from dG -> error
    with pytest.raises(ValueError):
        rf.getpairs(genes, snps, pd.DataFrame({"variant_id": ["sX"], "gene_id": ["A"]}))


def test_qvalue_bounded_and_ranks_with_probability():
    q = rf.qvalue([0.99, 0.5, 0.9, 0.1, 0.8])
    assert q.min() >= 0.0 and q.max() <= 1.0
    P = np.array([0.99, 0.5, 0.9, 0.1, 0.8])
    # the most probable edge must not have a larger q than the least probable one
    assert q[P.argmax()] <= q[P.argmin()]


def test_combine_formulas():
    res = {"p2": np.array([0.4]), "p3": np.array([0.5]),
           "p4": np.array([0.6]), "p5": np.array([0.8])}
    assert rf.combine(res, "orig")[0] == pytest.approx(0.5 * (0.4 * 0.8 + 0.6))
    assert rf.combine(res, "IV")[0] == pytest.approx(0.4 * 0.8)
    assert rf.combine(res, "mediation")[0] == pytest.approx(0.4 * 0.5)
    with pytest.raises(ValueError):
        rf.combine(res, "nope")


def test_dagfindr_breaks_cycle():
    # A->B (better q) is kept; B->A would close a cycle and is dropped
    dP = pd.DataFrame({"Source": ["A", "B"], "Target": ["B", "A"],
                       "Probability": [0.9, 0.8], "qvalue": [0.1, 0.2]})
    out, name2idx, G = rf.dagfindr_greedy_edges(dP)
    kept = {(r.Source, r.Target) for r in out[out.inDAG_greedy_edges].itertuples()}
    assert ("A", "B") in kept and ("B", "A") not in kept
    assert G.number_of_edges() == int(out.inDAG_greedy_edges.sum())
    assert list(out.columns) == ["Source", "Target", "Probability", "qvalue",
                                 "Source_idx", "Target_idx", "inDAG_greedy_edges"]


def test_write_all_edges_excludes_self_and_sorts(tmp_path):
    genes = ["g0", "g1", "g2", "g3"]
    colX = np.array([0, 1])                       # sources are g0, g1
    PP = np.array([[0.9, 0.8, 0.95, 0.2],
                   [0.1, 0.7, 0.85, 0.99]])
    self_mask = np.zeros(PP.shape, bool)
    self_mask[np.arange(2), colX] = True
    valid = ~self_mask
    qfull = np.zeros_like(PP)
    p = tmp_path / "T.grn_edges_all.csv.gz"
    n = rf.write_all_edges(str(p), PP, qfull, valid, genes, colX)
    d = pd.read_csv(p)
    assert n == 6 and len(d) == 6                 # 2 sources x 4 genes - 2 self = 6
    assert not (d.Source == d.Target).any()
    assert list(d.Probability) == sorted(d.Probability, reverse=True)


# --------------------------------------------------------------------------- end-to-end
def _write_inputs(d: Path, n=300, seed=0):
    """Two instrumented regulators: A->B and E->D are causal; C is null."""
    rng = np.random.default_rng(seed)
    snpA = rng.integers(0, 3, n).astype(float)
    snpE = rng.integers(0, 3, n).astype(float)
    A = 2.0 * snpA + rng.normal(0, 0.5, n)
    E = 2.0 * snpE + rng.normal(0, 0.5, n)
    B = 1.5 * A + rng.normal(0, 0.5, n)           # A drives B
    D = 1.5 * E + rng.normal(0, 0.5, n)           # E drives D
    C = rng.normal(0, 1.0, n)                      # null
    dX = pd.DataFrame({"A": A, "B": B, "C": C, "D": D, "E": E})
    dG = pd.DataFrame({"snpA": snpA.astype(int), "snpE": snpE.astype(int)})
    dE = pd.DataFrame({"variant_id": ["snpA", "snpE"], "gene_id": ["A", "E"]})
    dX.to_csv(d / "T.dX.csv", index=False)
    dG.to_csv(d / "T.dG.csv", index=False)
    dE.to_csv(d / "T.dE.csv", index=False)


def test_end_to_end_reconstruct(tmp_path):
    pytest.importorskip("findr")
    if not Path(rf.DEFAULT_LIBPATH).exists():
        pytest.skip("libfindr.so not available")
    _write_inputs(tmp_path)
    rf.reconstruct_grn("T", str(tmp_path), str(tmp_path), combination="orig", fdr=0.5,
                       method="kde", nth=2, na=2, libpath=rf.DEFAULT_LIBPATH, write_all=True)

    edges = pd.read_csv(tmp_path / "T.grn_edges.csv")
    allf = pd.read_csv(tmp_path / "T.grn_edges_all.csv.gz")
    qc = dict(r.split("\t") for r in (tmp_path / "T.grn_qc.tsv").read_text().splitlines()[1:])

    # schemas
    assert list(edges.columns) == ["Source", "Target", "Probability", "qvalue",
                                   "Source_idx", "Target_idx", "inDAG_greedy_edges"]
    assert list(allf.columns) == ["Source", "Target", "Probability", "qvalue"]
    # probabilities are valid and edges only start at instrumented sources
    assert edges.Probability.between(0, 1).all() and allf.Probability.between(0, 1).all()
    assert set(edges.Source) <= {"A", "E"}
    # master table: every non-self edge (2 sources x 5 genes - 2 self = 8), sorted
    assert len(allf) == 8
    assert list(allf.Probability) == sorted(allf.Probability, reverse=True)
    # findr recovers causality: the true causal edges outrank the null targets
    pa = dict(zip(zip(allf.Source, allf.Target), allf.Probability))
    assert pa[("A", "B")] > pa[("A", "C")]
    assert pa[("E", "D")] > pa[("E", "C")]
    # kept network contains the causal edges and is acyclic
    kept = set(zip(edges.Source, edges.Target))
    assert ("A", "B") in kept and ("E", "D") in kept
    import networkx as nx
    assert nx.is_directed_acyclic_graph(nx.DiGraph(list(kept)))
    # QC internally consistent
    assert int(qc["regulators"]) <= 2 and int(qc["edges"]) == int(edges.inDAG_greedy_edges.sum())
