#!/usr/bin/env python3
"""Combined validation figure, reading its inputs from results_findr_py/.

 (a) full per-tissue GRN overlaps with OmniPath over a degree-preserving null;
 (b) significant genes vs the GWAS sign-flip negative control;
 (c) the full FinnGen R11 TWAS replicates discovery -- left: discovery vs FinnGen
     z_TWAS over the discovery-significant genes, coronary heart disease (CHD) vs a
     height negative control; right: sign concordance by predictor class.

Panel (c) is the identical association test re-run in FinnGen over every tested
gene (analysis/validation/finngen_full_twas.py), not a lookup of the discovery hits.
It degrades gracefully (placeholder) if its inputs are not yet built.

Output: results_findr_py/figures/fig_flux_validation.{pdf,png}
"""
import os
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

RESULTS = Path(__file__).resolve().parents[2] / "results_findr_py"
GR = f"{RESULTS}/omnipath_validation"                     # omnipath_global_validation.csv
VD = f"{RESULTS}/omnipath_validation"
NEG = f"{VD}/negative_control_CAD_aragam2022.tsv"
FT_SUMMARY = f"{RESULTS}/finngen_twas/finngen_full_twas_summary.tsv"
FT_PERPAIR = f"{RESULTS}/finngen_twas/finngen_full_twas_perpair.tsv"
OUT = f"{RESULTS}/figures"
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
    ax.set_ylabel("OmniPath overlap\n(fold over degree null)", fontsize=9)
    ax.set_title("(a) Full GRN overlaps with OmniPath", fontsize=10, fontweight="bold", loc="left")
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


def panel_finn_scatter(ax):
    """Discovery vs FinnGen z_TWAS over the discovery-significant genes: CHD (real
    replication) vs the height negative control. One point per gene-tissue pair."""
    p = pd.read_csv(FT_PERPAIR, sep="\t")
    s = p[p.padj_disc < 0.05]
    chd = s[s.finn_trait == "CHD"]; hgt = s[s.finn_trait == "HEIGHT"]
    r_chd = np.corrcoef(chd.z_disc, chd.z_finn)[0, 1]
    r_hgt = np.corrcoef(hgt.z_disc, hgt.z_finn)[0, 1]
    lim = float(np.ceil(np.abs(np.r_[chd.z_disc, chd.z_finn, hgt.z_finn]).max()))
    ax.axhline(0, lw=0.6, color="#ddd", zorder=0); ax.axvline(0, lw=0.6, color="#ddd", zorder=0)
    ax.plot([-lim, lim], [-lim, lim], ls=":", lw=0.8, color="#888", zorder=1)
    ax.scatter(hgt.z_disc, hgt.z_finn, s=7, c=C_CTRL, alpha=0.45, lw=0, zorder=2)
    ax.scatter(chd.z_disc, chd.z_finn, s=7, c=C_CHD, alpha=0.55, lw=0, zorder=3)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal", "box")
    ax.set_xlabel("discovery $z_{\\mathrm{TWAS}}$", fontsize=9)
    ax.set_ylabel("FinnGen $z_{\\mathrm{TWAS}}$", fontsize=9)
    ax.text(0.04, 0.96, f"CHD  $r={r_chd:.2f}$", transform=ax.transAxes, ha="left", va="top",
            fontsize=8.5, fontweight="bold", color=C_CHD)
    ax.text(0.04, 0.87, f"height  $r={r_hgt:.2f}$", transform=ax.transAxes, ha="left", va="top",
            fontsize=8.5, color="#6b6f73")
    ax.set_title("(c) FinnGen replication: effect size", fontsize=10, fontweight="bold", loc="left")
    ax.legend(handles=[Patch(color=C_CHD, label="CHD"), Patch(color=C_CTRL, label="height")],
              frameon=False, fontsize=7.5, loc="lower right", handlelength=1.2, handletextpad=0.4)
    ax.tick_params(labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)


def panel_finn_bars(ax):
    """Sign concordance with discovery by predictor class: CHD vs height, over the
    discovery-significant stable genes."""
    s = pd.read_csv(FT_SUMMARY, sep="\t")

    def conc(tr, cls):
        row = s[(s.finn_trait == tr) & (s.subset == f"disc_sig_stable::{cls}")]
        return 100 * float(row.sign_conc.iloc[0]), float(row.sign_p.iloc[0])

    classes = ["cis_only", "cis_trans", "trans_only"]
    lab = {"cis_only": "cis-only", "cis_trans": "cis+trans", "trans_only": "trans-only"}
    x = np.arange(len(classes)); w = 0.4
    for i, cls in enumerate(classes):
        cc, cp = conc("CHD", cls); hc, _ = conc("HEIGHT", cls)
        ax.bar(x[i] - w / 2, cc, w, color=C_CHD)
        ax.bar(x[i] + w / 2, hc, w, color=C_CTRL)
        ax.text(x[i] - w / 2, cc + 1.5, f"{cc:.0f}\n{stars(cp)}", ha="center", va="bottom",
                fontsize=8, fontweight="bold", color=C_CHD)
        ax.text(x[i] + w / 2, hc + 1.5, f"{hc:.0f}", ha="center", va="bottom",
                fontsize=7.5, color="#666")
    ax.axhline(50, ls="--", lw=1, color="#444")
    ax.set_xticks(x); ax.set_xticklabels([lab[c] for c in classes], fontsize=8.5)
    ax.set_ylabel("sign concordance (%)", fontsize=9); ax.set_ylim(0, 112)
    ax.set_title("(d) FinnGen replication: concordance", fontsize=10, fontweight="bold", loc="left")
    ax.legend(handles=[Patch(color=C_CHD, label="CHD"), Patch(color=C_CTRL, label="height"),
                       Line2D([0], [0], ls="--", color="#444", label="chance (50%)")],
              frameon=False, fontsize=7.5, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.34),
              columnspacing=1.0, handlelength=1.2, handletextpad=0.4)
    ax.spines[["top", "right"]].set_visible(False)


def panel_placeholder(ax):
    ax.text(0.5, 0.5, "(c) FinnGen replication\n(pending)", ha="center", va="center",
            fontsize=10, color="#999", transform=ax.transAxes)
    ax.set_title("(c) Replicates in FinnGen CHD", fontsize=10, fontweight="bold", loc="left")
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)


def main():
    neg = pd.read_csv(NEG, sep="\t").set_index("metric")
    fig = plt.figure(figsize=(17.5, 3.9))
    gs = gridspec.GridSpec(1, 4, width_ratios=[1.6, 0.8, 1.3, 1.0], wspace=0.42, figure=fig)
    panel_enrich(fig.add_subplot(gs[0]))
    panel_collapse(fig.add_subplot(gs[1]), neg.loc["unique_sig_genes"],
                   "(b) Significant genes vs null")
    if os.path.exists(FT_PERPAIR) and os.path.exists(FT_SUMMARY):
        panel_finn_scatter(fig.add_subplot(gs[2]))
        panel_finn_bars(fig.add_subplot(gs[3]))
    else:
        panel_placeholder(fig.add_subplot(gs[2]))
    fig.text(0.5, -0.04, "significance:  * $p<0.05$   ** $p<0.01$   *** $p<0.001$   "
             "(a: vs degree null;  b: vs sign-flip null;  d: vs 50%)",
             ha="center", va="top", fontsize=8, color="#333")
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/fig_flux_validation.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_flux_validation.pdf")


if __name__ == "__main__":
    main()
