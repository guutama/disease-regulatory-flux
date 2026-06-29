"""Unit tests for bin/fit_expression_models.py (channel-aware Bayesian expression models).

These exercise the real NumPyro fit and ArviZ PSIS-LOO scoring, so the whole module is
skipped where NumPyro or ArviZ is unavailable (for example in the lightweight CI run). A
small tissue is simulated in tmp_path with three genes:
  * geneA -- driven by a cis SNP and also carrying a (noise) trans feature -> class "both",
  * geneB -- driven only by a trans feature -> class "trans_only",
  * geneC -- a cis SNP unrelated to its expression -> class "cis_only", not predictable.
We test the per-gene branching and config set, model selection and the predictability flag,
the weights/metrics/stats contracts (for the default bayes_ridge), and that the horseshoe and
bslmm priors also run and select sensibly.
"""
import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("numpyro")
pytest.importorskip("arviz")

_SPEC = importlib.util.spec_from_file_location(
    "fit_expression_models",
    Path(__file__).resolve().parents[1] / "bin" / "fit_expression_models.py")
fm = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fm)

N = 50
MCMC = ["--mcmc-warmup", "120", "--mcmc-samples", "120"]


def _write(path: Path, rows):
    path.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")


def build(d: Path):
    rng = np.random.default_rng(0)
    samples = [f"S{i}" for i in range(1, N + 1)]
    csA1, csA2, csC, trA, trB = (rng.integers(0, 3, size=N).astype(float) for _ in range(5))

    def std(x):
        return (x - x.mean()) / (x.std() + 1e-9)

    geneA = 1.8 * std(csA1) + 0.2 * rng.standard_normal(N)     # cis-driven
    geneB = 1.8 * std(trB) + 0.2 * rng.standard_normal(N)      # trans-driven
    geneC = rng.standard_normal(N)                             # noise

    _write(d / "expr.tsv", [["gene_id"] + samples,
                            ["geneA"] + list(np.round(geneA, 4)),
                            ["geneB"] + list(np.round(geneB, 4)),
                            ["geneC"] + list(np.round(geneC, 4))])
    _write(d / "geno.tsv", [["variant_id"] + samples,
                            ["csA1"] + list(csA1.astype(int)),
                            ["csA2"] + list(csA2.astype(int)),
                            ["csC"] + list(csC.astype(int)),
                            ["trA"] + list(trA.astype(int)),
                            ["trB"] + list(trB.astype(int))])
    _write(d / "cis.tsv", [["gene_id", "variant_id"],
                           ["geneA", "csA1"], ["geneA", "csA2"], ["geneC", "csC"]])
    _write(d / "trans.tsv", [["gene_id", "hop", "source_gene_id", "variant_id"],
                             ["geneA", "1", "geneX", "trA"],
                             ["geneB", "1", "geneR", "trB"]])


def run(d: Path, method="bayes_ridge"):
    build(d)
    fm.main(["--expression", str(d / "expr.tsv"), "--genotype", str(d / "geno.tsv"),
             "--cis-pruned", str(d / "cis.tsv"), "--trans-features", str(d / "trans.tsv"),
             "--tissue", "T", "--outdir", str(d), "--method", method] + MCMC)
    return {
        "metrics": d / "T.expr_model_metrics.tsv.gz",
        "weights": d / "T.expr_model_weights.tsv.gz",
        "selected": d / "T.expr_model_selected.tsv",
        "stats": d / "T.expr_model_stats.tsv",
    }


def read_gz(path: Path):
    return [ln.split("\t") for ln in gzip.decompress(path.read_bytes()).decode().splitlines()]


def read_tsv(path: Path):
    return [ln.split("\t") for ln in path.read_text().splitlines()]


def metrics_by_gene(path: Path):
    rows = read_gz(path)
    head = rows[0]
    return {r[0]: dict(zip(head, r)) for r in rows[1:]}


@pytest.fixture(scope="module")
def ridge(tmp_path_factory):
    d = tmp_path_factory.mktemp("ridge")
    return run(d, "bayes_ridge")


def test_gene_classes_and_configs(ridge):
    m = metrics_by_gene(ridge["metrics"])
    assert m["geneA"]["gene_class"] == "both"
    assert m["geneB"]["gene_class"] == "trans_only"
    assert m["geneC"]["gene_class"] == "cis_only"
    assert m["geneA"]["loo_r2_cis"] != "" and m["geneA"]["loo_r2_trans"] != "" \
        and m["geneA"]["loo_r2_cistrans"] != ""
    assert m["geneB"]["loo_r2_cis"] == "" and m["geneB"]["loo_r2_trans"] != ""
    assert m["geneC"]["loo_r2_trans"] == "" and m["geneC"]["loo_r2_cis"] != ""


def test_predictability(ridge):
    m = metrics_by_gene(ridge["metrics"])
    assert m["geneA"]["predictable"] == "True"        # cis signal
    assert m["geneB"]["predictable"] == "True"        # trans signal
    assert m["geneC"]["predictable"] == "False"       # noise
    assert m["geneB"]["best_config"] == "trans_only"


def test_weights_and_selected_contract(ridge):
    w = read_gz(ridge["weights"])
    assert w[0] == ["gene_id", "method", "config", "channel", "variant_id", "weight"]
    for r in w[1:]:
        assert r[1] == "bayes_ridge"
        float(r[5])
    csa = {r[4] for r in w[1:] if r[0] == "geneA" and r[2] == "cis_only"}
    assert csa == {"csA1", "csA2"}
    sel = read_tsv(ridge["selected"])
    assert sel[0] == ["gene_id", "gene_class", "best_config", "best_loo_r2",
                      "best_sigma_g", "best_intercept", "n_cis", "n_trans"]
    assert {r[0] for r in sel[1:]} == {"geneA", "geneB"}


def test_stats(ridge):
    stats = {r[0]: r[1] for r in read_tsv(ridge["stats"])[1:]}
    assert stats["method"] == "bayes_ridge"
    assert stats["n_genes"] == "3"
    assert stats["n_predictable"] == "2"
    assert stats["n_both"] == "1" and stats["n_trans_only"] == "1" and stats["n_cis_only"] == "1"


def test_horseshoe_runs(tmp_path):
    m = metrics_by_gene(run(tmp_path, "horseshoe")["metrics"])
    assert m["geneA"]["gene_class"] == "both"
    assert m["geneA"]["predictable"] == "True"
    assert m["geneC"]["predictable"] == "False"


def test_bslmm_runs(tmp_path):
    outs = run(tmp_path, "bslmm")
    m = metrics_by_gene(outs["metrics"])
    assert m["geneB"]["best_config"] == "trans_only"
    assert all(r[1] == "bslmm" for r in read_gz(outs["weights"])[1:])
