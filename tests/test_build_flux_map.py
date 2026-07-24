"""Unit tests for bin/build_flux_map.py (Step 10: disease-regulatory flux map).

The flux map attributes each disease gene's trans signal to its upstream regulators. Every
contributing ancestor is attributed DIRECTLY to the disease gene it feeds -- a hop-1 GRN parent
(hop=1) or a hop-2 grandparent (hop=2) -- so the per-regulator fluxes into a gene sum to
z_trans(g). An ancestor that delivers no weighted, GWAS-matched SNP is omitted. The output is one
CSV row per regulator->disease-gene attribution:
  regulator, regulator_gene_id, target, target_gene_id, tissue, hop, flux
We test:
  * a contributing parent edge (parent -> disease gene, hop 1),
  * a grandparent attributed directly to the disease gene (hop 2),
  * that a non-contributing ancestor is omitted,
  * the exact recovery sum_A flux(A -> g) = z_trans(g),
  * that only significant genes are decomposed.
"""
import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pytest

_SPEC = importlib.util.spec_from_file_location(
    "build_flux_map", Path(__file__).resolve().parents[1] / "bin" / "build_flux_map.py")
fm = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fm)

S = [f"S{i}" for i in range(8)]
D1 = [0, 1, 2, 0, 1, 2, 1, 0]      # s1 -- cis-eQTL of parent P1
D2 = [2, 2, 1, 0, 0, 1, 1, 2]      # s2 -- cis-eQTL of parent P2 (non-contributing: no weight)
D3 = [0, 0, 1, 1, 0, 1, 0, 2]      # s3 -- cis-eQTL of grandparent GP (reaches g via P1)
SD1, SD3 = np.array(D1, float).std(), np.array(D3, float).std()
SIGMA_G = 2.0
T1 = 0.5 * SD1 * 3.0 / SIGMA_G      # flux(P1 -> g)
T3 = -0.3 * SD3 * 2.0 / SIGMA_G     # flux(GP -> g)
Z_TRANS = T1 + T3                   # P2 contributes nothing

GENO = [["variant_id", "chromosome", "position", "ref", "alt"] + S,
        ["s1", "1", "101", "A", "G"] + D1,
        ["s2", "1", "102", "C", "T"] + D2,
        ["s3", "2", "201", "G", "A"] + D3]

# trans weights of geneG: s1 (from P1) and s3 (from GP) are weighted; s2 (P2) is not.
# Only the selected config (cis_trans) is used; the trans_only decoy must be ignored.
WEIGHTS = [["gene_id", "method", "config", "channel", "variant_id", "weight"],
           ["geneG", "horseshoe", "cis_trans",  "trans", "s1", "0.5"],
           ["geneG", "horseshoe", "cis_trans",  "trans", "s3", "-0.3"],
           ["geneG", "horseshoe", "trans_only", "trans", "s1", "99.0"]]   # decoy

SELECTED = [["gene_id", "gene_class", "best_config", "best_loo_r2", "best_elpd", "n_cis", "n_trans"],
            ["geneG", "both", "cis_trans", "0.1", "-1", "0", "2"]]

# trans features: g's parents P1 (via s1) and P2 (via s2); g's grandparent GP (via s3).
# P1's own parent is GP (so the GP -> P1 -> g path exists).
TRANS = [["gene_id", "hop", "source_gene_id", "variant_id"],
         ["geneG", "1", "P1", "s1"],
         ["geneG", "1", "P2", "s2"],
         ["geneG", "2", "GP", "s3"],
         ["P1",    "1", "GP", "s3"]]

ASSOC = [["gene_id", "config", "n_snp", "n_cis", "n_trans", "n_gwas", "z_twas", "z_cis",
          "z_trans", "sigma_g", "p_value", "p_adj", "tissue", "trait"],
         ["geneG", "cis_trans", "2", "0", "2", "2", f"{Z_TRANS:.10g}", "0",
          f"{Z_TRANS:.10g}", str(SIGMA_G), "1e-6", "0.01", "T", "CADtest"],
         ["geneN", "trans_only", "1", "0", "1", "1", "0.3", "0",
          "0.3", "1.0", "0.5", "0.4", "T", "CADtest"]]

GWAS = [["variant_id", "z", "beta", "se", "p_value", "n"],
        ["s1", "3.0", "0.3", "0.1", "1e-3", "1000"],
        ["s2", "1.0", "0.1", "0.1", "0.3", "1000"],
        ["s3", "2.0", "0.2", "0.1", "1e-2", "1000"]]


def _write(path: Path, rows):
    text = "\n".join("\t".join(map(str, r)) for r in rows) + "\n"
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def run(d: Path):
    _write(d / "geno.tsv.gz", GENO)
    _write(d / "w.tsv.gz", WEIGHTS)
    _write(d / "trans.tsv.gz", TRANS)
    _write(d / "assoc.tsv", ASSOC)
    _write(d / "gwas.tsv.gz", GWAS)
    _write(d / "sel.tsv", SELECTED)
    fm.main(["--association", str(d / "assoc.tsv"), "--weights", str(d / "w.tsv.gz"),
             "--selected", str(d / "sel.tsv"),
             "--trans-features", str(d / "trans.tsv.gz"), "--genotype", str(d / "geno.tsv.gz"),
             "--gwas", str(d / "gwas.tsv.gz"), "--tissue", "T", "--trait", "CADtest",
             "--outdir", str(d)])
    out = d / "flux_T_CADtest.csv"
    rows = [ln.split(",") for ln in out.read_text().splitlines()]
    header = rows[0]
    return [dict(zip(header, r)) for r in rows[1:]]


def _by_edge(rows):
    return {(r["regulator"], r["target"]): r for r in rows}


def test_parent_edge_flux(tmp_path):
    e = _by_edge(run(tmp_path))
    p1 = e[("P1", "geneG")]            # parent -> disease gene
    assert p1["hop"] == "1" and p1["target_gene_id"] == "geneG"
    assert float(p1["flux"]) == pytest.approx(T1)


def test_grandparent_attributed_directly(tmp_path):
    e = _by_edge(run(tmp_path))
    # the grandparent is attributed directly to the disease gene at hop 2
    gp = e[("GP", "geneG")]
    assert gp["hop"] == "2"
    assert float(gp["flux"]) == pytest.approx(T3)
    assert ("GP", "P1") not in e           # not routed through the intermediate parent


def test_noncontributing_ancestor_omitted(tmp_path):
    e = _by_edge(run(tmp_path))
    assert ("P2", "geneG") not in e        # P2 delivers no weighted SNP, so no edge


def test_flux_recovers_z_trans(tmp_path):
    rows = run(tmp_path)
    total = sum(float(r["flux"]) for r in rows if r["target"] == "geneG")
    assert total == pytest.approx(Z_TRANS)


def test_only_significant_genes_decomposed(tmp_path):
    targets = {r["target"] for r in run(tmp_path)}
    assert targets == {"geneG"}        # geneN (p_adj = 0.4) is not a disease gene
