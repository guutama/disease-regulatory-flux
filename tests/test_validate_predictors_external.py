"""Unit tests for the external predictor-validation script.

validate_predictors_external.py applies the published S9 weights to an external cohort's dosage
matrix to predict expression, then correlates against measured expression -- an out-of-sample
check of the discovery LOO-R2 with no retraining. These tests cover allele harmonisation, the
per-gene scoring, and the end-to-end file wiring.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pd = pytest.importorskip("pandas")

_SPEC = importlib.util.spec_from_file_location(
    "validate_predictors_external",
    Path(__file__).resolve().parents[1] / "analysis" / "validation"
    / "validate_predictors_external.py")
vpe = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(vpe)


def test_harmonise_flip_and_mismatch():
    alleles = {"v1": ("A", "G"), "v2": ("C", "T"), "v3": ("G", "A"), "v4": ("A", "T")}
    raw = {
        "v1": ("A", "G", np.array([0.0, 1.0, 2.0])),      # aligned -> as is
        "v2": ("T", "C", np.array([0.0, 1.0, 2.0])),      # swapped -> 2 - dosage
        "v3": (None, None, np.array([1.0, 1.0, 1.0])),    # no allele info -> assumed aligned
        "v4": ("A", "C", np.array([0.0, 1.0, 2.0])),      # different allele -> dropped
    }
    oriented, n_flip, n_mm = vpe.harmonise(raw, alleles)
    assert set(oriented) == {"v1", "v2", "v3"}            # v4 dropped
    assert np.allclose(oriented["v1"], [0.0, 1.0, 2.0])
    assert np.allclose(oriented["v2"], [2.0, 1.0, 0.0])   # flipped
    assert n_flip == 1 and n_mm == 1


def test_impute_fills_and_drops():
    mat = np.array([[0.0, np.nan, 2.0],      # one gap -> filled with row mean (1.0)
                    [np.nan, np.nan, np.nan]])  # all missing -> dropped
    out, keep = vpe._impute(mat)
    assert keep.tolist() == [True, False]
    assert np.allclose(out, [[0.0, 1.0, 2.0]])


def _weights():
    return pd.DataFrame([
        # g1: two cis SNPs, both present externally
        {"gene": "GENE1", "gene_id": "ENSG1.1", "tissue": "AOR", "config": "cis_only",
         "rs_id": "rs1", "weight": 1.0},
        {"gene": "GENE1", "gene_id": "ENSG1.1", "tissue": "AOR", "config": "cis_only",
         "rs_id": "rs2", "weight": 2.0},
        # g2: its only SNP is absent from the external dosage -> unscored
        {"gene": "GENE2", "gene_id": "ENSG2.1", "tissue": "AOR", "config": "trans_only",
         "rs_id": "rs3", "weight": 1.0},
    ])


def test_compute_scores_present_gene_and_reports_missing():
    dosage = pd.DataFrame(
        {"s1": [0.0, 1.0], "s2": [1.0, 1.0], "s3": [2.0, 0.0], "s4": [1.0, 2.0], "s5": [0.0, 1.0]},
        index=["rs1", "rs2"])
    yhat = 1.0 * dosage.loc["rs1"] + 2.0 * dosage.loc["rs2"]      # [2,3,2,5,2]
    expr = pd.DataFrame([yhat.to_numpy()], index=["ENSG1"], columns=dosage.columns)
    loo = {"ENSG1.1": (0.5, True), "ENSG2.1": (0.2, True)}

    out = vpe.compute(_weights(), dosage, expr, loo, min_snp=1)
    g1 = out[out.gene_id == "ENSG1.1"].iloc[0]
    assert g1["n_snp_used"] == 2 and g1["n_snp_missing"] == 0 and g1["n_samples"] == 5
    assert g1["r2_external"] == pytest_approx(1.0)               # predictor reproduces expression
    assert g1["loo_r2"] == 0.5 and bool(g1["predictable"])
    # g2's SNP is missing externally and it is not measured -> no score, missing counted
    g2 = out[out.gene_id == "ENSG2.1"].iloc[0]
    assert g2["n_snp_used"] == 0 and g2["n_snp_missing"] == 1
    assert np.isnan(g2["r2_external"])


def pytest_approx(x, tol=1e-9):
    class _A:
        def __eq__(self, other):
            return abs(other - x) <= tol
    return _A()


def test_main_end_to_end(tmp_path):
    dos = tmp_path / "dosage.tsv"
    dos.write_text("rs_id\ts1\ts2\ts3\ts4\ts5\n"
                   "rs1\t0\t1\t2\t1\t0\n"
                   "rs2\t1\t1\t0\t2\t1\n")
    expr = tmp_path / "expr.tsv"
    expr.write_text("gene_id\ts1\ts2\ts3\ts4\ts5\n"
                    "ENSG1.1\t2\t3\t2\t5\t2\n")
    s9 = tmp_path / "s9.csv"
    pd.DataFrame([
        {"gene": "GENE1", "gene_id": "ENSG1.1", "tissue": "AOR", "config": "cis_only",
         "rs_id": "rs1", "ref": "A", "alt": "G", "weight": 1.0},
        {"gene": "GENE1", "gene_id": "ENSG1.1", "tissue": "AOR", "config": "cis_only",
         "rs_id": "rs2", "ref": "C", "alt": "T", "weight": 2.0},
    ]).to_csv(s9, index=False)
    s1 = tmp_path / "s1.csv"
    pd.DataFrame([{"gene_id": "ENSG1.1", "tissue": "AOR", "best_loo_r2": 0.4,
                   "predictable": "True"}]).to_csv(s1, index=False)
    out = tmp_path / "out.csv"

    vpe.main(["--weights", str(s9), "--loo", str(s1), "--dosage", str(dos),
              "--expression", str(expr), "--tissue", "AOR", "--out", str(out)])
    res = pd.read_csv(out)
    row = res[res.gene_id == "ENSG1.1"].iloc[0]
    assert row["n_snp_used"] == 2
    assert abs(row["r2_external"] - 1.0) < 1e-9
    assert row["loo_r2"] == 0.4
