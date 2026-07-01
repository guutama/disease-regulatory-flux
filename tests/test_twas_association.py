"""Unit tests for bin/twas_association.py (Stage 9: summary-statistic TWAS association).

For each gene's selected predictor we combine its ALT-dosage weights with the ALT-aligned
GWAS z to form z_TWAS = sum_j w_j * sd_j * z_gwas_j / sigma_g, where sd_j is the SNP's
genotype standard deviation and sigma_g is the predicted-expression spread. The statistic
splits additively into a cis and a trans component. We test:
  * the numerator / sigma_g / z_TWAS arithmetic against an independent NumPy computation,
  * the additive cis + trans channel split,
  * the two-sided p-value and per-tissue BH-FDR,
  * that a gene with no GWAS-matched SNP (or a degenerate predictor) is dropped.
"""
import gzip
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

_SPEC = importlib.util.spec_from_file_location(
    "twas_association", Path(__file__).resolve().parents[1] / "bin" / "twas_association.py")
ta = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ta)

S = [f"S{i}" for i in range(8)]

# genotype_012: variant_id, chromosome, position, ref, alt, then samples
V1 = [0, 1, 2, 0, 1, 2, 1, 0]      # cis SNP of geneA, also geneB
V2 = [2, 2, 1, 0, 0, 1, 1, 2]      # trans SNP of geneA
V3 = [0, 0, 1, 1, 0, 1, 0, 0]      # cis SNP of geneC (no GWAS match -> geneC dropped)
GENO = [["variant_id", "chromosome", "position", "ref", "alt"] + S,
        ["v1", "1", "101", "A", "G"] + V1,
        ["v2", "2", "201", "C", "T"] + V2,
        ["v3", "3", "301", "G", "A"] + V3]

# weights: the file carries ALL fitted configs per gene; only the selected one is used.
# The 99.0 decoy rows are in non-selected configs and must be ignored.
WEIGHTS = [["gene_id", "method", "config", "channel", "variant_id", "weight"],
           ["geneA", "horseshoe", "cis_trans",  "cis",   "v1", "0.5"],
           ["geneA", "horseshoe", "cis_trans",  "trans", "v2", "-0.3"],
           ["geneA", "horseshoe", "cis_only",   "cis",   "v1", "99.0"],   # decoy
           ["geneA", "horseshoe", "trans_only", "trans", "v2", "99.0"],   # decoy
           ["geneB", "horseshoe", "cis_only",   "cis",   "v1", "0.8"],
           ["geneB", "horseshoe", "cis_trans",  "cis",   "v1", "99.0"],   # decoy
           ["geneC", "horseshoe", "cis_only",   "cis",   "v3", "0.4"],
           ["geneD", "horseshoe", "cis_only",   "cis",   "v1", "0.7"]]    # not predictable

# selected config per predictable gene (geneD is absent -> not predictable -> not tested)
SELECTED = [["gene_id", "gene_class", "best_config", "best_loo_r2", "best_elpd", "n_cis", "n_trans"],
            ["geneA", "both", "cis_trans", "0.1", "-1", "1", "1"],
            ["geneB", "cis", "cis_only",  "0.1", "-1", "1", "0"],
            ["geneC", "cis", "cis_only",  "0.1", "-1", "1", "0"]]

# harmonised GWAS: variant_id, z, beta, se, p_value, n  (v3 absent -> geneC has no match)
GWAS = [["variant_id", "z", "beta", "se", "p_value", "n"],
        ["v1", "3.0", "0.30", "0.10", "1e-3", "1000"],
        ["v2", "-1.5", "-0.15", "0.10", "0.13", "1000"]]


def _write(path: Path, rows):
    text = "\n".join("\t".join(map(str, r)) for r in rows) + "\n"
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


def run(d: Path):
    _write(d / "geno.tsv.gz", GENO)
    _write(d / "w.tsv.gz", WEIGHTS)
    _write(d / "gwas.tsv.gz", GWAS)
    _write(d / "sel.tsv", SELECTED)
    ta.main(["--weights", str(d / "w.tsv.gz"), "--genotype", str(d / "geno.tsv.gz"),
             "--gwas", str(d / "gwas.tsv.gz"), "--selected", str(d / "sel.tsv"),
             "--tissue", "T", "--trait", "CADtest", "--outdir", str(d)])
    out = d / "association_T_CADtest.tsv"
    rows = [ln.split("\t") for ln in out.read_text().splitlines()]
    header = rows[0]
    return [dict(zip(header, r)) for r in rows[1:]]


def _expected(weights):
    """Independent NumPy computation of z_twas / z_cis / z_trans for one gene."""
    G = {"v1": np.array(V1, float), "v2": np.array(V2, float), "v3": np.array(V3, float)}
    z = {"v1": 3.0, "v2": -1.5}
    sd = {v: G[v].std() for v in G}
    pred = sum(w * G[v] for _, v, w in weights)
    sigma_g = pred.std()
    num = {"cis": 0.0, "trans": 0.0}
    for ch, v, w in weights:
        if v in z:
            num[ch] += w * sd[v] * z[v]
    z_cis, z_trans = num["cis"] / sigma_g, num["trans"] / sigma_g
    return z_cis, z_trans, z_cis + z_trans, sigma_g


def test_ztwas_and_channel_split(tmp_path):
    res = {r["gene_id"]: r for r in run(tmp_path)}
    z_cis, z_trans, z_twas, sigma_g = _expected(
        [("cis", "v1", 0.5), ("trans", "v2", -0.3)])
    a = res["geneA"]
    assert a["config"] == "cis_trans"
    assert int(a["n_cis"]) == 1 and int(a["n_trans"]) == 1 and int(a["n_gwas"]) == 2
    assert float(a["sigma_g"]) == pytest.approx(sigma_g)
    assert float(a["z_cis"]) == pytest.approx(z_cis)
    assert float(a["z_trans"]) == pytest.approx(z_trans)
    assert float(a["z_twas"]) == pytest.approx(z_twas)
    # additive split
    assert float(a["z_cis"]) + float(a["z_trans"]) == pytest.approx(float(a["z_twas"]))


def test_cis_only_has_zero_trans(tmp_path):
    res = {r["gene_id"]: r for r in run(tmp_path)}
    b = res["geneB"]
    assert b["config"] == "cis_only"
    assert float(b["z_trans"]) == pytest.approx(0.0)
    assert float(b["z_cis"]) == pytest.approx(float(b["z_twas"]))


def test_pvalue_two_sided(tmp_path):
    res = {r["gene_id"]: r for r in run(tmp_path)}
    a = res["geneA"]
    z = float(a["z_twas"])
    assert float(a["p_value"]) == pytest.approx(math.erfc(abs(z) / math.sqrt(2)))


def test_unmatched_gene_dropped(tmp_path):
    res = {r["gene_id"]: r for r in run(tmp_path)}
    # geneC's only SNP (v3) has no GWAS record -> no association emitted
    assert "geneC" not in res


def test_bh_fdr_monotone(tmp_path):
    res = run(tmp_path)
    pairs = sorted((float(r["p_value"]), float(r["p_adj"])) for r in res)
    # BH-adjusted p-values are non-decreasing in raw p and never below the raw p
    adj = [a for _, a in pairs]
    assert all(adj[i] <= adj[i + 1] + 1e-12 for i in range(len(adj) - 1))
    assert all(a >= p - 1e-12 for p, a in pairs)


def test_uses_only_selected_config_and_predictable(tmp_path):
    res = {r["gene_id"]: r for r in run(tmp_path)}
    assert "geneD" not in res                       # absent from selected -> not tested
    assert res["geneA"]["config"] == "cis_trans"    # the selected config is reported
    assert res["geneB"]["config"] == "cis_only"
