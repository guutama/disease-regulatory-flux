#!/usr/bin/env python3
"""Combined validation figure (4 panels):
 (a) flux-map edges reproduce OmniPath interactions over a degree-preserving null
     (per tissue; directed and undirected);
 (b) significant genes collapse under the GWAS sign-flip negative control;
 (c) OmniPath-confirmed disease edges collapse under the same control;
 (d) the predicted direction replicates in an independent GWAS (FinnGen R11 CHD)
     but not in a height negative control, by predictor class.
Panels (b,c) show observed (red) vs the permutation null (grey: mean bar, whisker
to max). Significance is shown with stars throughout (key at the foot).

Inputs (per-edge / per-gene SUMMARY tables only -- no individual-level data), under RESULTS:
  <GR>/<tissue>.grn_edges.tsv                          flux-map edges vs OmniPath (per tissue)
  <VD>/negative_control_<trait>.tsv                    GWAS sign-flip negative control
  <VD>/finngen_replication_summary.tsv                 independent-GWAS (FinnGen CHD) replication

Output: fig_flux_validation.{pdf,png} in the shared figures directory.

Run:  python analysis/figures/plot_flux_validation.py

Input locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS); figures
are written to results/figures/ by default (override with FLUX_FIGURES). Style: Arial, 300 DPI.
"""
import os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS = Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_cv"))
GR = f"{RESULTS}/horseshoe_alltargets/flux/flux_grn"
VD = f"{RESULTS}/horseshoe_alltargets/omnipath_validation"
NEG = f"{VD}/negative_control_CAD_aragam2022.tsv"
FINN = f"{VD}/finngen_replication_summary.tsv"
# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
OUT = str(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))
os.makedirs(OUT, exist_ok=True)
TIS = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
C_OBS, C_NULL, C_CHD, C_CTRL = "#c1121f", "#9aa0a6", "#c1121f", "#9aa0a6"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 9, "pdf.fonttype": 42, "savefig.dpi": 300})


def stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def panel_enrich(ax):
    d = pd.read_csv(f"{GR}/omnipath_global_validation.csv")
    d = d[d.dataset == "GRN"]                                  # the FULL per-tissue GRN
    dt = d[d.network.isin(TIS)].set_index("network").reindex(TIS).reset_index()   # all 7 tissues
    x = np.arange(len(dt)); w = 0.38
    for key, name, col, off in [("directed", "directed", "#C44E52", -0.5),
                                ("undirected", "undirected", "#4C72B0", 0.5)]:
        ax.bar(x + off * w, dt[f"{key}_enrich"], w, label=name, color=col, edgecolor="white", lw=0.5)
        for xi, v, p in zip(x + off * w, dt[f"{key}_enrich"], dt[f"{key}_p"]):
            if p < 0.05:
                ax.text(xi, v + 0.05, stars(p), ha="center", va="bottom", fontsize=9)
    ax.axhline(1.0, ls="--", lw=0.9, color="#888", zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(dt.network, fontsize=8)
    ax.set_ylim(0, max(dt.directed_enrich.max(), dt.undirected_enrich.max()) * 1.22)
    ax.set_ylabel("OmniPath reproduction\n(fold over degree null)", fontsize=9)
    ax.set_title("(a) Full GRN reproduces OmniPath", fontsize=10, fontweight="bold", loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)


def panel_collapse(ax, row, title):
    ax.bar(["sign-flip\nnull"], [row.null_mean], width=0.55, color=C_NULL, edgecolor="white", zorder=2)
    ax.errorbar(["sign-flip\nnull"], [row.null_mean], yerr=[[0], [row.null_max - row.null_mean]],
                fmt="none", ecolor="#555", capsize=5, lw=1.2, zorder=3)
    ax.bar(["observed"], [row.observed], width=0.55, color=C_OBS, edgecolor="white", zorder=2)
    top = max(row.observed, row.null_max)
    ax.set_ylim(0, top * 1.22)
    ax.text(1, row.observed + top * 0.03, f"{int(row.observed):,}\n{stars(row.emp_p)}", ha="center",
            va="bottom", fontsize=8.5, fontweight="bold", color=C_OBS)
    ax.text(0, row.null_max + top * 0.03, f"max {int(row.null_max)}", ha="center",
            va="bottom", fontsize=7.5, color="#555")
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    ax.tick_params(labelsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)


def panel_replication(ax):
    rep = pd.read_csv(FINN, sep="\t")
    chd = rep[rep.finn == "CHD"].set_index("class"); hgt = rep[rep.finn == "HEIGHT"].set_index("class")
    classes = ["cis_only", "cis_trans", "trans_only"]
    lab = {"cis_only": "cis-only", "cis_trans": "cis+trans", "trans_only": "trans-only"}
    x = np.arange(len(classes)); w = 0.4
    ax.bar(x - w / 2, [100 * chd.loc[c, "concordance"] for c in classes], w, color=C_CHD)
    ax.bar(x + w / 2, [100 * hgt.loc[c, "concordance"] for c in classes], w, color=C_CTRL)
    for xi, c in zip(x, classes):
        v = 100 * chd.loc[c, "concordance"]
        ax.text(xi - w / 2, v + 1.5, f"{v:.0f}\n{stars(chd.loc[c,'p'])}", ha="center",
                va="bottom", fontsize=8, fontweight="bold", color=C_CHD)
        ax.text(xi + w / 2, 100 * hgt.loc[c, "concordance"] + 1.5,
                f"{100*hgt.loc[c,'concordance']:.0f}", ha="center", va="bottom", fontsize=7.5, color="#666")
    ax.axhline(50, ls="--", lw=1, color="#444")
    ax.set_xticks(x); ax.set_xticklabels([lab[c] for c in classes], fontsize=8.5)
    ax.set_ylabel("sign concordance (%)", fontsize=9); ax.set_ylim(0, 112)
    ax.set_title("(d) Replicates in FinnGen CHD", fontsize=10, fontweight="bold", loc="left")
    ax.legend(handles=[Patch(color=C_CHD, label="CHD"), Patch(color=C_CTRL, label="height"),
                       Line2D([0], [0], ls="--", color="#444", label="chance (50%)")],
              frameon=False, fontsize=7.5, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.34),
              columnspacing=1.0, handlelength=1.2, handletextpad=0.4)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    neg = pd.read_csv(NEG, sep="\t").set_index("metric")
    fig = plt.figure(figsize=(16.5, 3.8))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1.7, 0.85, 0.85, 1.3], wspace=0.34, figure=fig)
    panel_enrich(fig.add_subplot(gs[0]))
    panel_collapse(fig.add_subplot(gs[1]), neg.loc["unique_sig_genes"],
                   "(b) Genes vs null")
    panel_collapse(fig.add_subplot(gs[2]), neg.loc["omnipath_confirmed_directed_edges"],
                   "(c) Edges vs null")
    panel_replication(fig.add_subplot(gs[3]))
    fig.text(0.5, -0.04, "significance:  * $p<0.05$   ** $p<0.01$   *** $p<0.001$   "
             "(a: vs degree null;  b, c: vs sign-flip null;  d: vs 50%)",
             ha="center", va="top", fontsize=8, color="#333")
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/fig_flux_validation.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_flux_validation.pdf")


if __name__ == "__main__":
    main()
