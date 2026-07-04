"""Unit tests for analysis/validation/twas_from_gwas.py (external summary-statistic TWAS).

These check the core algebra of `compute` -- the cis/trans decomposition, the sigma_g
normalisation, and the skipping of SNPs that lack an sd or a GWAS z -- plus an end-to-end
main() run that writes a table and joins the discovery comparison.
"""
import importlib.util
import math
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

_SPEC = importlib.util.spec_from_file_location(
    "twas_from_gwas",
    Path(__file__).resolve().parents[1] / "analysis" / "validation" / "twas_from_gwas.py")
tw = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(tw)


def _weights():
    # G1 (cis_trans): rs1 cis, rs2 trans, rs3 trans-with-blank-sd (skipped).
    # G2 (cis_only): rs4 cis, but the GWAS has no z for it -> G2 dropped.
    return pd.DataFrame([
        {"gene": "GENE1", "gene_id": "G1", "tissue": "T1", "config": "cis_trans", "channel": "cis",   "rs_id": "rs1", "weight": 2.0, "sd": 1.0},
        {"gene": "GENE1", "gene_id": "G1", "tissue": "T1", "config": "cis_trans", "channel": "trans", "rs_id": "rs2", "weight": 3.0, "sd": 0.5},
        {"gene": "GENE1", "gene_id": "G1", "tissue": "T1", "config": "cis_trans", "channel": "trans", "rs_id": "rs3", "weight": 9.0, "sd": float("nan")},
        {"gene": "GENE2", "gene_id": "G2", "tissue": "T1", "config": "cis_only",  "channel": "cis",   "rs_id": "rs4", "weight": 1.0, "sd": 1.0},
    ])


def test_compute_decomposition_and_normalisation():
    gwas_z = {"rs1": 1.0, "rs2": 2.0, "rs3": 5.0}     # rs4 absent -> G2 dropped
    sigma_g = {("G1", "T1"): 2.0}
    out = tw.compute(_weights(), gwas_z, sigma_g)

    assert list(out["gene_id"]) == ["G1"]             # G2 has no usable GWAS SNP
    r = out.iloc[0]
    # z_cis = 2*1*1 / 2 = 1.0 ; z_trans = 3*0.5*2 / 2 = 1.5 ; rs3 skipped (blank sd)
    assert r["z_cis"] == pytest.approx(1.0)
    assert r["z_trans"] == pytest.approx(1.5)
    assert r["z_twas"] == pytest.approx(2.5)
    assert (r["n_snp"], r["n_cis"], r["n_trans"], r["n_gwas"]) == (2, 1, 1, 2)
    assert r["p"] == pytest.approx(math.erfc(2.5 / math.sqrt(2)))


def test_compute_skips_zero_and_missing_sigma_g():
    out = tw.compute(_weights(), {"rs1": 1.0, "rs2": 2.0}, {("G1", "T1"): 0.0})
    assert out.empty                                   # sigma_g == 0 -> gene dropped


def test_main_end_to_end_with_discovery(tmp_path):
    w = tmp_path / "s9.csv"
    _weights().to_csv(w, index=False)
    s2 = tmp_path / "s2.csv"
    pd.DataFrame([{"gene_id": "G1", "tissue": "T1", "sigma_g": 2.0, "z_twas": 4.0}]).to_csv(s2, index=False)
    gwas = tmp_path / "gwas.tsv"
    pd.DataFrame([{"rs_id": "rs1", "z": 1.0}, {"rs_id": "rs2", "z": 2.0}]).to_csv(gwas, sep="\t", index=False)
    out = tmp_path / "out.csv"

    tw.main(["--weights", str(w), "--sigma-g", str(s2), "--gwas", str(gwas),
             "--out", str(out), "--discovery", str(s2)])

    got = pd.read_csv(out)
    assert list(got.columns)[:14] == tw.OUT_COLS
    row = got.iloc[0]
    assert row["z_twas"] == pytest.approx(2.5)
    # discovery z_twas = 4.0 (same sign) -> concordant
    assert row["z_discovery"] == 4.0 and bool(row["concordant"]) is True
    assert 0.0 <= row["p_adj"] <= 1.0
