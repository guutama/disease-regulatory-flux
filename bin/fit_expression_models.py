#!/usr/bin/env python3
"""
Fit a channel-aware Bayesian expression model for every gene in one tissue.

Each gene's expression is predicted from two genetic channels: its own cis-eQTL SNPs (the
cis channel) and the cis-eQTL SNPs of its upstream regulators in the network (the trans
channel). For each gene the model is fit in up to three configurations and the best one is
kept:

  cis_only    the cis channel alone,
  trans_only  the trans channel alone,
  cis_trans   both channels together.

A gene with both channels is fit in all three; a gene with only one channel is fit in that
one. The configuration with the highest model evidence is selected, and the gene counts as
predictable when the selected model's leave-one-out R^2 reaches a threshold.

Method (--method):
  bayes_ridge  (default) channel-aware empirical-Bayes ridge: a separate Gaussian prior
               variance per channel, fit by EM, with an exact closed-form leave-one-out R^2
               and log-evidence. Fast --- no sampling and no train/test split.
  horseshoe    a regularised-horseshoe model fit by MCMC (heavier; selected by ELPD). Set up
               by --method horseshoe; requires NumPyro.

Inputs (one tissue):
  --expression     harmonised gene-by-sample expression matrix (the response per gene);
  --genotype       0/1/2 hard-call genotype matrix (variant_id + sample columns);
  --cis-pruned     LD-pruned per-gene cis-SNP set (gene_id, variant_id);
  --trans-features per-gene trans-feature SNPs (gene_id, hop, source_gene_id, variant_id).

Output (in --outdir), for tissue <tissue>:
  <tissue>.expr_model_metrics.tsv.gz  one row per gene: per-configuration LOO-R^2, evidence
                                      and intercept, the gene class, the selected config and
                                      whether the gene is predictable;
  <tissue>.expr_model_weights.tsv.gz  long posterior-mean raw-dosage weights:
                                      gene_id, method, config, channel, variant_id, weight;
  <tissue>.expr_model_selected.tsv    the predictable genes and their selected config;
  <tissue>.expr_model_stats.tsv       tissue-level counts and mean LOO-R^2 per channel.

Column names and the field delimiters are options, so a new dataset is supported by adjusting
configuration rather than editing code.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict
from typing import TextIO

import numpy as np

PRED_THRESHOLD = 0.01
LOG2PI = float(np.log(2.0 * np.pi))


# ----------------------------------------------------------------------------- io helpers
def _open(path: str, mode: str = "rt") -> TextIO:
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _split(line: str, delim: str) -> list[str]:
    return line.rstrip("\n").split(delim)


def _col(header: list[str], name: str, where: str) -> int:
    try:
        return header.index(name)
    except ValueError:
        raise SystemExit(f"The {where} file is missing the required column '{name}'.")


# -------------------------------------------------------------------------------- loaders
def load_gene_variants(path: str, delim: str, gene_col: str, variant_col: str):
    """gene_id -> list of variant ids (first-seen order, de-duplicated)."""
    out: dict[str, list[str]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        gi = _col(header, gene_col, "feature")
        vi = _col(header, variant_col, "feature")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            key = (f[gi], f[vi])
            if key not in seen:
                seen.add(key)
                out[f[gi]].append(f[vi])
    return out


def load_expression(path: str, delim: str, gene_col: str, genes: set[str]):
    """Return (samples, gene_id -> expression vector) for the requested genes."""
    expr: dict[str, np.ndarray] = {}
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        if header[0] != gene_col:
            raise SystemExit(f"The first column of the expression matrix must be "
                             f"'{gene_col}', but it is '{header[0]}'.")
        samples = header[1:]
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            if f[0] in genes:
                expr[f[0]] = np.array(f[1:], dtype=np.float64)
    return samples, expr


def load_genotype(path: str, delim: str, variant_col: str, wanted: set[str]):
    """Return (samples, variant_id -> dosage vector) for the requested variants."""
    geno: dict[str, np.ndarray] = {}
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        if header[0] != variant_col:
            raise SystemExit(f"The first column of the genotype matrix must be "
                             f"'{variant_col}', but it is '{header[0]}'.")
        samples = header[1:]
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            if f[0] in wanted and f[0] not in geno:
                geno[f[0]] = np.array(f[1:], dtype=np.float64)
    return samples, geno


# --------------------------------------------------------------- empirical-Bayes ridge fit
def eb_ridge_loo(Xstd: np.ndarray, channel: np.ndarray, y: np.ndarray, em_iters: int = 15):
    """Channel-aware empirical-Bayes ridge on standardised X and centred y.

    Each coefficient has a Gaussian prior with a per-channel precision (lambda_cis,
    lambda_trans), fit by EM together with the noise variance. Returns the posterior-mean
    coefficients (in standardised-X space), the exact leave-one-out R^2 (from the ridge
    leverage), the in-sample R^2 and the log model evidence."""
    n, p = Xstd.shape
    Xty = Xstd.T @ y
    XtX = Xstd.T @ Xstd
    masks = (channel == 0, channel == 1)
    lam = np.ones(2)
    eye = np.eye(p)
    for _ in range(em_iters):
        Lam = lam[channel]
        Ainv = np.linalg.inv(XtX + np.diag(Lam) + 1e-9 * eye)
        beta = Ainv @ Xty
        diag_invA = np.diag(Ainv)
        gamma = np.clip(1.0 - Lam * diag_invA, 0.0, 1.0)
        for c, mask in enumerate(masks):
            if not mask.any():
                continue
            num = float(gamma[mask].sum())
            den = float(np.sum(beta[mask] ** 2))
            lam[c] = num / den if (den > 0 and num > 0) else max(lam[c], 1e-8)

    Lam = lam[channel]
    A = XtX + np.diag(Lam) + 1e-9 * eye
    Ainv = np.linalg.inv(A)
    beta = Ainv @ Xty
    yhat = Xstd @ beta
    h = np.einsum("ij,ij->i", Xstd @ Ainv, Xstd)        # leverage H_ii
    denom = np.clip(1.0 - h, 1e-8, None)
    sst = float(np.sum((y - y.mean()) ** 2))
    loo_resid = (y - yhat) / denom
    loo_r2 = float(1 - np.sum(loo_resid ** 2) / sst) if sst > 0 else float("nan")
    in_r2 = float(1 - np.sum((y - yhat) ** 2) / sst) if sst > 0 else float("nan")

    gamma_sum = float(np.sum(np.clip(1.0 - Lam * np.diag(Ainv), 0.0, 1.0)))
    sigma2 = max(float(np.sum((y - yhat) ** 2)) / max(n - gamma_sum, 1.0), 1e-8)
    sign, logdetA = np.linalg.slogdet(A)
    logdetB = float(logdetA) - float(np.sum(np.log(Lam)))
    yBy = float(y @ y) - float(Xty @ beta)
    log_evidence = -0.5 * (n * LOG2PI + n * np.log(sigma2) + logdetB + yBy / sigma2)
    return dict(beta_std=beta, loo_r2=loo_r2, in_r2=in_r2, log_evidence=log_evidence)


def standardise(cols: list[np.ndarray]):
    """Stack genotype columns, drop constant ones, return (Xstd, mean, sd, keep_mask)."""
    M = np.column_stack(cols).astype(np.float64)
    mu = M.mean(0)
    sd = M.std(0)
    keep = sd > 0
    return (M[:, keep] - mu[keep]) / sd[keep], mu[keep], sd[keep], keep


def topk_by_corr(cols: list[np.ndarray], y: np.ndarray, k: int) -> list[int]:
    """Indices of the k columns most correlated (absolute) with y."""
    M = np.column_stack(cols).astype(np.float64)
    Mc = M - M.mean(0)
    yc = y - y.mean()
    denom = np.sqrt((Mc ** 2).sum(0) * (yc @ yc)) + 1e-12
    r = np.abs((Mc.T @ yc) / denom)
    return sorted(range(M.shape[1]), key=lambda j: -r[j])[:k]


# --------------------------------------------------------------- regularised horseshoe fit
def fit_horseshoe(Xstd: np.ndarray, channel: np.ndarray, y: np.ndarray,
                  *, warmup: int, samples: int, seed: int,
                  slab_scale: float = 2.0, slab_df: float = 4.0, tau0: float = 0.1):
    """Channel-aware regularised-horseshoe model fit by MCMC (NumPyro).

    Each coefficient has a horseshoe prior with a per-channel global scale; the slab
    regularises the largest coefficients (Piironen & Vehtari 2017). Returns the posterior-mean
    coefficients (standardised-X space), the in-sample R^2, an approximate leave-one-out R^2
    (truncated importance sampling over the posterior draws) and a WAIC elpd used to select
    among configurations (the role the ridge log-evidence plays for bayes_ridge)."""
    import jax
    import jax.numpy as jnp
    import numpyro
    import numpyro.distributions as dist
    from numpyro.infer import MCMC, NUTS

    n, p = Xstd.shape
    n_ch = int(channel.max()) + 1

    def model(X, ch, obs=None):
        sigma = numpyro.sample("sigma", dist.HalfNormal(1.0))
        tau = numpyro.sample("tau", dist.HalfCauchy(tau0 * jnp.ones(n_ch)))
        lam = numpyro.sample("lam", dist.HalfCauchy(jnp.ones(p)))
        c2 = numpyro.sample("c2", dist.InverseGamma(slab_df / 2., (slab_df / 2.) * slab_scale ** 2))
        z = numpyro.sample("z", dist.Normal(0., 1.).expand([p]).to_event(1))
        tau_c = tau[ch]
        lam_t = jnp.sqrt(c2 * lam ** 2 / (c2 + tau_c ** 2 * lam ** 2))
        beta = numpyro.deterministic("beta", z * tau_c * lam_t)
        numpyro.sample("obs", dist.Normal(X @ beta, sigma), obs=obs)

    numpyro.set_host_device_count(1)
    mcmc = MCMC(NUTS(model), num_warmup=warmup, num_samples=samples,
                num_chains=1, progress_bar=False)
    mcmc.run(jax.random.PRNGKey(seed), jnp.asarray(Xstd), jnp.asarray(channel),
             obs=jnp.asarray(y))
    post = mcmc.get_samples()
    beta_s = np.asarray(post["beta"])            # (S, p)
    sigma_s = np.asarray(post["sigma"])          # (S,)
    beta_mean = beta_s.mean(0)
    mu_s = beta_s @ Xstd.T                        # (S, n) predicted means per draw
    yhat = Xstd @ beta_mean
    sst = float(np.sum((y - y.mean()) ** 2))
    in_r2 = float(1 - np.sum((y - yhat) ** 2) / sst) if sst > 0 else float("nan")

    # pointwise log-likelihood L (S, n); WAIC elpd for model selection
    resid = y[None, :] - mu_s
    L = -0.5 * (LOG2PI + 2 * np.log(sigma_s)[:, None] + (resid ** 2) / (sigma_s ** 2)[:, None])
    s = L.shape[0]
    lppd = np.log(np.mean(np.exp(L - L.max(0)), axis=0)) + L.max(0)
    elpd_waic = float(np.sum(lppd - np.var(L, axis=0)))

    # leave-one-out predictive mean by truncated importance sampling (weights ~ 1/likelihood)
    logr = -L
    logr -= logr.max(0)
    r = np.exp(logr)
    w = np.minimum(r, r.mean(0) * np.sqrt(s))
    yhat_loo = (w * mu_s).sum(0) / np.clip(w.sum(0), 1e-12, None)
    loo_r2 = float(1 - np.sum((y - yhat_loo) ** 2) / sst) if sst > 0 else float("nan")
    return dict(beta_std=beta_mean, loo_r2=loo_r2, in_r2=in_r2, log_evidence=elpd_waic)


# ----------------------------------------------------------------------------------- main
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--expression", required=True, help="input: harmonised expression matrix")
    p.add_argument("--genotype", required=True, help="input: 0/1/2 hard-call genotype matrix")
    p.add_argument("--cis-pruned", required=True, help="input: LD-pruned per-gene cis SNPs")
    p.add_argument("--trans-features", required=True, help="input: per-gene trans-feature SNPs")
    p.add_argument("--tissue", required=True, help="tissue label, used to name the outputs")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--method", choices=["bayes_ridge", "horseshoe"], default="bayes_ridge",
                   help="expression model (default: bayes_ridge)")
    p.add_argument("--min-cis", type=int, default=1,
                   help="minimum cis SNPs for a gene to get a cis channel (default: 1)")
    p.add_argument("--max-trans", type=int, default=200,
                   help="cap on trans features per gene, kept by correlation (default: 200)")
    p.add_argument("--pred-threshold", type=float, default=PRED_THRESHOLD,
                   help="selected LOO-R^2 at/above which a gene is predictable (default: 0.01)")
    p.add_argument("--mcmc-warmup", type=int, default=500,
                   help="horseshoe: MCMC warm-up iterations (default: 500)")
    p.add_argument("--mcmc-samples", type=int, default=500,
                   help="horseshoe: MCMC sampling iterations (default: 500)")
    p.add_argument("--mcmc-seed", type=int, default=0,
                   help="horseshoe: MCMC random seed (default: 0)")
    # format options
    p.add_argument("--delimiter", default="\t", help="delimiter of the tab inputs (default: tab)")
    p.add_argument("--gene-col", default="gene_id", help="gene-id column (default: gene_id)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column (default: variant_id)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def gene_class(n_cis: int, n_trans: int) -> str:
    if n_cis > 0 and n_trans > 0:
        return "both"
    return "cis_only" if n_cis > 0 else "trans_only"


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    cis_map = load_gene_variants(args.cis_pruned, d, args.gene_col, args.variant_col)
    trans_map = load_gene_variants(args.trans_features, d, args.gene_col, args.variant_col)
    cis_cand = {g for g, v in cis_map.items() if len(v) >= args.min_cis}
    trans_cand = {g for g, v in trans_map.items() if len(v) >= 1}
    genes = sorted(cis_cand | trans_cand)
    if not genes:
        raise SystemExit("No gene has a cis or trans channel, so nothing can be fit.")

    universe: set[str] = set()
    for g in genes:
        universe.update(cis_map.get(g, ()))
        universe.update(trans_map.get(g, ()))

    expr_samples, expr = load_expression(args.expression, d, args.gene_col, set(genes))
    geno_samples, geno = load_genotype(args.genotype, d, args.variant_col, universe)

    geno_pos = {s: i for i, s in enumerate(geno_samples)}
    common = [s for s in expr_samples if s in geno_pos]
    if not common:
        raise SystemExit("The expression and genotype matrices share no samples.")
    e_pos = {s: i for i, s in enumerate(expr_samples)}
    e_idx = np.array([e_pos[s] for s in common])
    g_idx = np.array([geno_pos[s] for s in common])

    os.makedirs(args.outdir, exist_ok=True)
    weights_path = os.path.join(args.outdir, f"{args.tissue}.expr_model_weights.tsv.gz")
    metrics_path = os.path.join(args.outdir, f"{args.tissue}.expr_model_metrics.tsv.gz")
    selected_path = os.path.join(args.outdir, f"{args.tissue}.expr_model_selected.tsv")

    wf = _open(weights_path, "wt")
    wf.write("gene_id\tmethod\tconfig\tchannel\tvariant_id\tweight\n")
    metrics: list[dict] = []

    def fit_config(gene: str, config: str, variants: list[str], channels: list[str],
                   yc: np.ndarray, ymu: float):
        cols = [geno[v][g_idx] for v in variants]
        Xstd, mu, sd, keep = standardise(cols)
        if Xstd.shape[1] == 0:
            return None
        kept = [i for i, k in enumerate(keep) if k]
        ch = np.array([0 if channels[i] == "cis" else 1 for i in kept])
        if args.method == "horseshoe":
            r = fit_horseshoe(Xstd, ch, yc, warmup=args.mcmc_warmup,
                              samples=args.mcmc_samples, seed=args.mcmc_seed)
        else:
            r = eb_ridge_loo(Xstd, ch, yc)
        w_raw = r["beta_std"] / sd
        b0 = ymu - float(np.sum(w_raw * mu))
        for i, w in zip(kept, w_raw):
            wf.write(f"{gene}\t{args.method}\t{config}\t{channels[i]}\t{variants[i]}\t{w}\n")
        n_cis = sum(1 for i in kept if channels[i] == "cis")
        return dict(gene_id=gene, config=config, n_cis=n_cis, n_trans=len(kept) - n_cis,
                    loo_r2=r["loo_r2"], in_r2=r["in_r2"], log_evidence=r["log_evidence"],
                    intercept=b0)

    for gene in genes:
        if gene not in expr:
            continue
        y = expr[gene][e_idx]
        ymu = float(y.mean())
        yc = y - ymu
        cis_v = [v for v in cis_map.get(gene, ()) if v in geno]
        cis_set = set(cis_v)
        trans_v = [v for v in trans_map.get(gene, ()) if v in geno and v not in cis_set]
        if len(trans_v) > args.max_trans:
            cols = [geno[v][g_idx] for v in trans_v]
            trans_v = [trans_v[i] for i in topk_by_corr(cols, yc, args.max_trans)]
        has_cis = len(cis_v) >= args.min_cis

        if has_cis:
            m = fit_config(gene, "cis_only", cis_v, ["cis"] * len(cis_v), yc, ymu)
            if m:
                metrics.append(m)
            if trans_v:
                m = fit_config(gene, "trans_only", trans_v, ["trans"] * len(trans_v), yc, ymu)
                if m:
                    metrics.append(m)
                allv = cis_v + trans_v
                ch = ["cis"] * len(cis_v) + ["trans"] * len(trans_v)
                m = fit_config(gene, "cis_trans", allv, ch, yc, ymu)
                if m:
                    metrics.append(m)
        elif trans_v:
            m = fit_config(gene, "trans_only", trans_v, ["trans"] * len(trans_v), yc, ymu)
            if m:
                metrics.append(m)
    wf.close()

    write_outputs(metrics, args, metrics_path, selected_path)


def write_outputs(metrics: list[dict], args, metrics_path: str, selected_path: str):
    """Pivot the per-gene-config rows to one row per gene, select the best config by evidence,
    flag predictable genes, and write the metrics, selected and stats files."""
    by_gene: dict[str, dict[str, dict]] = defaultdict(dict)
    for m in metrics:
        by_gene[m["gene_id"]][m["config"]] = m
    short = {"cis_only": "cis", "trans_only": "trans", "cis_trans": "cistrans"}

    metric_cols = (["gene_id", "gene_class", "n_cis", "n_trans"]
                   + [f"loo_r2_{c}" for c in short.values()]
                   + [f"logE_{c}" for c in short.values()]
                   + ["best_config", "best_loo_r2", "predictable"])
    rows = []
    for gene in sorted(by_gene):
        cfgs = by_gene[gene]
        n_cis = max((m["n_cis"] for m in cfgs.values()), default=0)
        n_trans = max((m["n_trans"] for m in cfgs.values()), default=0)
        row = {"gene_id": gene, "gene_class": gene_class(n_cis, n_trans),
               "n_cis": n_cis, "n_trans": n_trans}
        for cfg, s in short.items():
            row[f"loo_r2_{s}"] = cfgs[cfg]["loo_r2"] if cfg in cfgs else ""
            row[f"logE_{s}"] = cfgs[cfg]["log_evidence"] if cfg in cfgs else ""
        best = max(cfgs.values(), key=lambda m: m["log_evidence"])
        row["best_config"] = best["config"]
        row["best_loo_r2"] = best["loo_r2"]
        row["predictable"] = best["loo_r2"] >= args.pred_threshold
        rows.append(row)

    with _open(metrics_path, "wt") as fh:
        fh.write("\t".join(metric_cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in metric_cols) + "\n")

    with _open(selected_path, "wt") as fh:
        fh.write("gene_id\tgene_class\tbest_config\tbest_loo_r2\tn_cis\tn_trans\n")
        for row in rows:
            if row["predictable"]:
                fh.write(f"{row['gene_id']}\t{row['gene_class']}\t{row['best_config']}\t"
                         f"{row['best_loo_r2']}\t{row['n_cis']}\t{row['n_trans']}\n")

    write_stats(rows, args)


def write_stats(rows: list[dict], args):
    stats_path = os.path.join(args.outdir, f"{args.tissue}.expr_model_stats.tsv")
    classes = ("cis_only", "trans_only", "both")
    configs = ("cis_only", "trans_only", "cis_trans")
    pred = [r for r in rows if r["predictable"]]

    def mean_loo(key: str) -> float:
        vals = [float(r[key]) for r in rows if r[key] != ""]
        return float(np.mean(vals)) if vals else float("nan")

    lines = [("tissue", args.tissue), ("method", args.method), ("n_genes", len(rows))]
    for cls in classes:
        lines.append((f"n_{cls}", sum(1 for r in rows if r["gene_class"] == cls)))
    lines.append(("n_predictable", len(pred)))
    for cls in classes:
        lines.append((f"n_predictable_{cls}", sum(1 for r in pred if r["gene_class"] == cls)))
    for cfg in configs:
        lines.append((f"n_selected_{cfg}", sum(1 for r in pred if r["best_config"] == cfg)))
    for s in ("cis", "trans", "cistrans"):
        lines.append((f"loo_r2_{s}_mean", mean_loo(f"loo_r2_{s}")))
    with _open(stats_path, "wt") as fh:
        fh.write("metric\tvalue\n")
        for k, v in lines:
            fh.write(f"{k}\t{v}\n")
    sys.stderr.write(f"[fit_expression_models] tissue={args.tissue} method={args.method} "
                     f"genes={len(rows)} predictable={len(pred)}\n")


if __name__ == "__main__":
    main()
