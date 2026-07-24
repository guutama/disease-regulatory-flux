#!/usr/bin/env python3
"""Pairwise agreement of the three shrinkage priors on per-gene predictive accuracy.

For every gene-model predictable by at least one prior (selected-configuration
LOO-R^2 >= 0.01), the held-out R^2 of the three priors is compared pairwise with a
seaborn pairplot. Points hug the diagonal, i.e. the priors rank and score genes
almost identically.

Output: results_findr_py/figures/fig_prior_pairplot.{pdf,png}
"""
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "results_findr_py"
OUT = RF / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PRIORS = [("Bayesian ridge", "expression_models_ridge"),
          ("Horseshoe", "expression_models"),
          ("BSLMM", "expression_models_bslmm")]


def load(subdir, r2col):
    fs = sorted(glob.glob(f"{RF/subdir}/*.expr_model_metrics.tsv.gz"))
    d = pd.concat([pd.read_csv(f, sep="\t", usecols=["gene_id", "best_loo_r2"])
                   .assign(tissue=Path(f).name.split(".")[0]) for f in fs], ignore_index=True)
    return d.rename(columns={"best_loo_r2": r2col})


def main():
    cols = [p[0] for p in PRIORS]
    m = None
    for name, sub in PRIORS:
        d = load(sub, name)
        m = d if m is None else m.merge(d, on=["tissue", "gene_id"], how="inner")
    m = m.loc[m[cols].max(axis=1) >= 0.01, cols]      # predictable by >=1 prior
    print(f"{len(m):,} gene-models predictable by >=1 prior")

    sns.set_theme(style="ticks", font="DejaVu Sans")
    g = sns.pairplot(m, vars=cols, diag_kind="hist", corner=False,
                     plot_kws=dict(s=6, alpha=0.15, edgecolor="none", color="#4C72B0"),
                     diag_kws=dict(color="#4C72B0", bins=60))
    lim = (min(0.0, m.min().min()), np.ceil(m.max().max() * 10) / 10)
    for i in range(len(cols)):
        for j in range(len(cols)):
            ax = g.axes[i, j]
            if i != j:
                ax.plot(lim, lim, ls=":", lw=0.9, color="#c1121f", zorder=5)  # y = x
                r = m[cols[j]].corr(m[cols[i]])
                ax.text(0.05, 0.92, f"r={r:.4f}", transform=ax.transAxes, fontsize=8,
                        color="#c1121f", fontweight="bold")
                ax.set_xlim(lim); ax.set_ylim(lim)
    g.figure.suptitle("Per-gene held-out $R^2$: pairwise agreement of the three priors", y=1.02,
                      fontsize=11, fontweight="bold")
    for ext in ("pdf", "png"):
        g.figure.savefig(OUT / f"fig_prior_pairplot.{ext}", dpi=300, bbox_inches="tight")
    plt.close(g.figure)
    print(f"wrote {OUT}/fig_prior_pairplot.pdf")


if __name__ == "__main__":
    main()
