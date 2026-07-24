#!/usr/bin/env python3
"""Supplementary benchmark of the three channel-aware shrinkage priors on the findr harness.

Two panels: (a) per-channel held-out accuracy (mean LOO-R^2 for cis, trans, cis+trans),
showing the priors are near-identical; (b) training cost per job (wall time and peak RAM),
showing where they differ. Accuracy is computed from the expression-model metrics; the cost
numbers are the measured per-job mean wall time and peak resident set (SLURM sacct).

Output: results_findr_py/figures/fig_prior_benchmark.{pdf,png}
"""
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RF = ROOT / "results_findr_py"
OUT = RF / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PRIORS = [("Bayesian ridge", "expression_models_ridge", "#4C72B0"),
          ("Horseshoe",      "expression_models",       "#C44E52"),
          ("BSLMM",          "expression_models_bslmm", "#55A868")]
# measured training cost per job (SLURM sacct): mean wall (min), peak RAM (GB)
COST = {"Bayesian ridge": (11.5, 2.4), "Horseshoe": (27.3, 10.5), "BSLMM": (27.2, 2.6)}
CHANNELS = [("loo_r2_cis", "cis"), ("loo_r2_trans", "trans"), ("loo_r2_cistrans", "cis+trans")]

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 9, "pdf.fonttype": 42, "savefig.dpi": 300})


def per_channel(subdir):
    d = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(glob.glob(f"{RF/subdir}/*.expr_model_metrics.tsv.gz"))],
                  ignore_index=True)
    return [pd.to_numeric(d[c], errors="coerce").dropna().mean() for c, _ in CHANNELS]


def main():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.5, 3.6), gridspec_kw={"width_ratios": [1.35, 1.0]})

    # (a) accuracy per channel
    x = np.arange(len(CHANNELS)); w = 0.26
    for i, (name, sub, col) in enumerate(PRIORS):
        axA.bar(x + (i - 1) * w, per_channel(sub), w, label=name, color=col, edgecolor="white", lw=0.4)
    axA.set_xticks(x); axA.set_xticklabels([c for _, c in CHANNELS])
    axA.set_ylabel("mean held-out $R^2$")
    axA.set_title("(a) Accuracy is near-identical", fontsize=10, fontweight="bold", loc="left")
    axA.legend(frameon=False, fontsize=8, loc="upper left")
    axA.spines[["top", "right"]].set_visible(False)

    # (b) cost: wall time (bars) + peak RAM (points on twin axis)
    names = [p[0] for p in PRIORS]; cols = [p[2] for p in PRIORS]
    xi = np.arange(len(names))
    axB.bar(xi, [COST[n][0] for n in names], 0.55, color=cols, edgecolor="white", lw=0.4)
    axB.set_ylabel("wall time / job (min)")
    axB.set_xticks(xi); axB.set_xticklabels([n.replace("Bayesian ", "Bayesian\n") for n in names], fontsize=8)
    axB.set_title("(b) Cost differs", fontsize=10, fontweight="bold", loc="left")
    axB.spines[["top"]].set_visible(False)
    axR = axB.twinx()
    axR.plot(xi, [COST[n][1] for n in names], "o-", color="#333", lw=1.2, ms=6)
    axR.set_ylabel("peak RAM (GB)")
    axR.set_ylim(0, 12)
    for i, n in enumerate(names):
        axR.annotate(f"{COST[n][1]:.1f} GB", (xi[i], COST[n][1]), textcoords="offset points",
                     xytext=(0, 7), ha="center", fontsize=7.5, color="#333")
    axR.spines[["top"]].set_visible(False)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig_prior_benchmark.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_prior_benchmark.pdf")


if __name__ == "__main__":
    main()
