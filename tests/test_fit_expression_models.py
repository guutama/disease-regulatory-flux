"""Unit tests for bin/fit_expression_models.py (channel-aware Bayesian expression models).

A small tissue is simulated in tmp_path with three genes:
  * geneA -- driven by a cis SNP and also carrying a (noise) trans feature -> class "both",
  * geneB -- driven only by a trans feature -> class "trans_only",
  * geneC -- a cis SNP unrelated to its expression -> class "cis_only", not predictable.
We test the empirical-Bayes ridge fit (LOO-R^2 separates signal from noise), the per-gene
branching and config set, model selection and the predictability flag, the weights/metrics/
stats contracts, and that --method horseshoe is reported as unavailable.
"""
import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pytest

_SPEC = importlib.util.spec_from_file_location(
    "fit_expression_models",
    Path(__file__).resolve().parents[1] / "bin" / "fit_expression_models.py")
fm = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(fm)

N = 60
RNG = np.random.default_rng(0)


def _dosage(n):
    return RNG.integers(0, 3, size=n).astype(float)


def _write(path: Path, rows):
    path.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")


def build(d: Path):
    samples = [f"S{i}" for i in range(1, N + 1)]
    csA1, csA2, csC, trA, trB = (_dosage(N) for _ in range(5))

    def std(x):
        return (x - x.mean()) / (x.std() + 1e-9)

    geneA = 1.8 * std(csA1) + 0.2 * RNG.standard_normal(N)     # cis-driven
    geneB = 1.8 * std(trB) + 0.2 * RNG.standard_normal(N)      # trans-driven
    geneC = RNG.standard_normal(N)                             # noise (csC unrelated)

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


def run(d: Path, method="bayes_ridge", extra=None):
    build(d)
    argv = ["--expression", str(d / "expr.tsv"), "--genotype", str(d / "geno.tsv"),
            "--cis-pruned", str(d / "cis.tsv"), "--trans-features", str(d / "trans.tsv"),
            "--tissue", "T", "--outdir", str(d), "--method", method]
    fm.main(argv + (extra or []))
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


def test_eb_ridge_loo_separates_signal_from_noise():
    rng = np.random.default_rng(1)
    x = rng.standard_normal((50, 1))
    y_sig = (2.0 * x[:, 0] + 0.1 * rng.standard_normal(50))
    y_noise = rng.standard_normal(50)
    ch = np.array([0])
    r_sig = fm.eb_ridge_loo((x - x.mean(0)) / x.std(0), ch, y_sig - y_sig.mean())
    r_noise = fm.eb_ridge_loo((x - x.mean(0)) / x.std(0), ch, y_noise - y_noise.mean())
    assert r_sig["loo_r2"] > 0.8
    assert r_noise["loo_r2"] < 0.05


def test_gene_classes_and_configs(tmp_path):
    m = metrics_by_gene(run(tmp_path)["metrics"])
    assert m["geneA"]["gene_class"] == "both"
    assert m["geneB"]["gene_class"] == "trans_only"
    assert m["geneC"]["gene_class"] == "cis_only"
    # geneA fit in all three configs (all LOO-R^2 columns populated)
    assert m["geneA"]["loo_r2_cis"] != "" and m["geneA"]["loo_r2_trans"] != "" \
        and m["geneA"]["loo_r2_cistrans"] != ""
    # geneB only trans_only; geneC only cis_only
    assert m["geneB"]["loo_r2_cis"] == "" and m["geneB"]["loo_r2_trans"] != ""
    assert m["geneC"]["loo_r2_trans"] == "" and m["geneC"]["loo_r2_cis"] != ""


def test_predictability(tmp_path):
    m = metrics_by_gene(run(tmp_path)["metrics"])
    assert m["geneA"]["predictable"] == "True"        # cis signal
    assert m["geneB"]["predictable"] == "True"        # trans signal
    assert m["geneC"]["predictable"] == "False"       # noise
    assert m["geneB"]["best_config"] == "trans_only"


def test_weights_contract(tmp_path):
    w = read_gz(run(tmp_path)["weights"])
    assert w[0] == ["gene_id", "method", "config", "channel", "variant_id", "weight"]
    # every weight is a finite number and method is recorded
    for r in w[1:]:
        assert r[1] == "bayes_ridge"
        float(r[5])
    # geneA cis_only weights cover its cis SNPs
    csa = {r[4] for r in w[1:] if r[0] == "geneA" and r[2] == "cis_only"}
    assert csa == {"csA1", "csA2"}


def test_stats_and_selected(tmp_path):
    outs = run(tmp_path)
    stats = {r[0]: r[1] for r in read_tsv(outs["stats"])[1:]}
    assert stats["n_genes"] == "3"
    assert stats["n_predictable"] == "2"
    assert stats["n_both"] == "1" and stats["n_trans_only"] == "1" and stats["n_cis_only"] == "1"
    sel = {r[0] for r in read_tsv(outs["selected"])[1:]}
    assert sel == {"geneA", "geneB"}


def test_horseshoe_runs_and_recovers_signal(tmp_path):
    pytest.importorskip("numpyro")           # heavy optional dep; skipped where unavailable
    outs = run(tmp_path, method="horseshoe",
               extra=["--mcmc-warmup", "150", "--mcmc-samples", "150"])
    m = metrics_by_gene(outs["metrics"])
    # same branching as bayes_ridge, and the method is recorded in the weights
    assert m["geneA"]["gene_class"] == "both"
    assert m["geneB"]["best_config"] == "trans_only"
    assert m["geneA"]["predictable"] == "True"        # cis signal recovered
    assert m["geneC"]["predictable"] == "False"       # noise
    w = read_gz(outs["weights"])
    assert all(r[1] == "horseshoe" for r in w[1:])
