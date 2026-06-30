"""Unit tests for bin/build_flux_map.py (Step 10: disease-regulatory flux map).

The flux map is a path-preserving network. Its nodes are the disease genes (significant in
the association) plus their hop-1 GRN parents and hop-2 grandparents; its edges are the true
single-hop GRN edges along the paths down to a disease gene: parent -> disease_gene (hop 1)
and grandparent -> intermediate_parent (hop 2). Each edge carries the signed flux its source
regulator delivers to the disease gene, flux(A -> g); a grandparent reaching g through several
parents splits evenly across those edges, so the fluxes still sum to z_trans(g). Every parent
and grandparent appears even when it delivers no flux (blank `flux`). We test:
  * a contributing parent edge (parent -> disease gene),
  * a grandparent edge whose target is the intermediate parent (the child gene is added),
  * that a non-contributing ancestor is kept with a blank flux,
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
WEIGHTS = [["gene_id", "method", "config", "channel", "variant_id", "weight"],
           ["geneG", "horseshoe", "cis_trans", "trans", "s1", "0.5"],
           ["geneG", "horseshoe", "cis_trans", "trans", "s3", "-0.3"]]

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
    fm.main(["--association", str(d / "assoc.tsv"), "--weights", str(d / "w.tsv.gz"),
             "--trans-features", str(d / "trans.tsv.gz"), "--genotype", str(d / "geno.tsv.gz"),
             "--gwas", str(d / "gwas.tsv.gz"), "--tissue", "T", "--trait", "CADtest",
             "--outdir", str(d)])
    out = d / "flux_T_CADtest.tsv"
    rows = [ln.split("\t") for ln in out.read_text().splitlines()]
    header = rows[0]
    return [dict(zip(header, r)) for r in rows[1:]]


def _by_edge(rows):
    return {(r["source_gene"], r["target_gene"]): r for r in rows}


def test_parent_edge_flux(tmp_path):
    e = _by_edge(run(tmp_path))
    p1 = e[("P1", "geneG")]            # parent -> disease gene
    assert p1["hop"] == "1" and p1["disease_gene"] == "geneG"
    assert float(p1["flux"]) == pytest.approx(T1)


def test_grandparent_edge_targets_intermediate(tmp_path):
    e = _by_edge(run(tmp_path))
    # the grandparent edge points at the intermediate parent P1, not at geneG
    gp = e[("GP", "P1")]
    assert gp["hop"] == "2" and gp["disease_gene"] == "geneG"
    assert float(gp["flux"]) == pytest.approx(T3)
    assert ("GP", "geneG") not in e        # no collapsed grandparent -> disease shortcut


def test_noncontributing_ancestor_blank_flux(tmp_path):
    e = _by_edge(run(tmp_path))
    p2 = e[("P2", "geneG")]            # P2 is a parent but delivers no weighted SNP
    assert p2["flux"] == ""            # kept as a structural edge, flux blank


def test_flux_recovers_z_trans(tmp_path):
    rows = run(tmp_path)
    total = sum(float(r["flux"]) for r in rows
                if r["disease_gene"] == "geneG" and r["flux"] != "")
    assert total == pytest.approx(Z_TRANS)


def test_only_significant_genes_decomposed(tmp_path):
    diseases = {r["disease_gene"] for r in run(tmp_path)}
    assert diseases == {"geneG"}       # geneN (p_adj = 0.4) is not a disease gene
