#!/usr/bin/env python3
"""OmniPath-corroboration figures for the disease flux maps.

(1) fig_flux_omnipath_null: per-tissue and pooled fold-enrichment of OmniPath
    reproduction (directed and undirected) over a degree-preserving null.
    Reference line at 1x; * marks p < 0.05.

(2) fig_flux_omnipath_network (the supplementary figure): exemplar panels showing the
    grandparent->parent->disease cascade, with OmniPath-corroborated relations highlighted
    (red = OmniPath relation, blue = GRN relation) and, where the corroboration is indirect,
    the intermediate gene it relays through.

Inputs are the flux-network OmniPath overlap tables under RESULTS (flux/flux_grn/). Outputs:
fig_flux_omnipath_{null,network}.{pdf,png} in the shared figures directory.

Input locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS); figures
go to results/figures/ by default (override with FLUX_FIGURES). Style: Arial, 300 DPI.
"""
import os
from pathlib import Path

import numpy as np, pandas as pd, networkx as nx
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D

# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS = Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_cv"))
GR = f"{RESULTS}/horseshoe_alltargets/flux/flux_grn"
# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
OUT = str(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))
TIS = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
C_HUB, C_PARENT, C_GP = "#6a51a3", "#2171b5", "#9ecae1"
C_DIR, C_VIA, C_SHORT, C_STRUCT = "#c1121f", "#f08c00", "#1b9e77", "#c2c5cc"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "pdf.fonttype": 42, "savefig.dpi": 300})


def fig_null():
    d = pd.read_csv(f"{GR}/omnipath_global_validation.csv")
    d = d[d.dataset == "FLUX"].copy()
    # drop tissues with no exact OmniPath match at all (empty bars, e.g. Blood, LIV)
    d = d[(d.directed_enrich > 0) | (d.undirected_enrich > 0)]
    d = d.set_index("network").reindex([n for n in TIS if n in set(d.network)]).reset_index()
    lab = list(d.network)
    x = np.arange(len(d)); w = 0.38
    fig, ax = plt.subplots(figsize=(8.0, 3.4))
    for j, (key, name, col) in enumerate([("directed", "directed", "#C44E52"),
                                          ("undirected", "undirected", "#4C72B0")]):
        vals = d[f"{key}_enrich"].values
        ps = d[f"{key}_p"].values
        ax.bar(x + (j - 0.5) * w, vals, w, label=name, color=col, edgecolor="white", lw=0.5)
        for xi, v, p in zip(x + (j - 0.5) * w, vals, ps):
            if p < 0.05:
                ax.text(xi, v + 0.06, "*", ha="center", va="bottom", fontsize=12, color="#222")
    ax.axhline(1.0, ls="--", lw=0.9, color="#888", zorder=0)
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=0)
    ax.set_ylabel("OmniPath reproduction\n(fold over degree-preserving null)")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/fig_flux_omnipath_null.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_flux_omnipath_null.pdf")


def _layout(parents, gps_xmap, gps):
    Y_HUB, Y_P, Y_GP = 0.0, 1.7, 3.4                 # more vertical space between layers
    pos = {"__core__": np.array([0.0, Y_HUB])}
    sp = 2.3                                          # horizontal spacing per node
    n = len(parents)
    xs = [0.0] if n <= 1 else list(np.linspace(-sp * (n - 1) / 2, sp * (n - 1) / 2, n))
    px = {}
    for p, xx in zip(parents, xs):
        px[p] = float(xx); pos[p] = np.array([float(xx), Y_P])
    if gps:
        gtx = {g: (np.mean([px[p] for p in gps_xmap[g] if p in px]) if any(p in px for p in gps_xmap[g]) else 0.0)
               for g in gps}
        ordered = sorted(gps, key=lambda g: gtx[g])
        m = len(ordered)
        xs_g = [gtx[ordered[0]]] if m == 1 else list(np.linspace(-sp * (m - 1) / 2, sp * (m - 1) / 2, m))
        for g, xe in zip(ordered, xs_g):
            pos[g] = np.array([0.5 * gtx[g] + 0.5 * xe, Y_GP])
    return pos, px


def _draw_case(ax, mode, reg, core, inter):
    """One exemplar subplot: disease core at the bottom, regulator on top, optional
    grey intermediate. Red = OmniPath relation, blue = flux-map (our GRN) relation.
    mode: 'both' (both direct), 'opvia' (our direct / OmniPath via intermediate),
          'ourvia' (OmniPath direct / our 2-hop via intermediate)."""
    from matplotlib.patches import Circle
    C_REG, C_INT, C_CORE, C_OP, C_GRN = "#1f6fb2", "#bdbdbd", "#6a51a3", "#c1121f", "#2c6fbb"
    if mode == "both":
        P_REG, P_CORE = (0.0, 0.9), (0.0, -0.9)
        ax.add_patch(FancyArrowPatch(P_REG, P_CORE, arrowstyle="-|>", mutation_scale=14,
                     color=C_OP, lw=2.2, connectionstyle="arc3,rad=0.32", shrinkA=18, shrinkB=20, zorder=3))
        ax.add_patch(FancyArrowPatch(P_REG, P_CORE, arrowstyle="-|>", mutation_scale=12,
                     color=C_GRN, lw=1.8, connectionstyle="arc3,rad=-0.32", shrinkA=16, shrinkB=18, zorder=3))
    else:
        P_REG, P_INT, P_CORE = (-0.40, 0.85), (0.62, 0.0), (-0.40, -0.85)
        direct_col, via_col = (C_OP, C_GRN) if mode == "ourvia" else (C_GRN, C_OP)
        ax.add_patch(FancyArrowPatch(P_REG, P_CORE, arrowstyle="-|>", mutation_scale=15,
                     color=direct_col, lw=2.3, connectionstyle="arc3,rad=0.42", shrinkA=18, shrinkB=20, zorder=3))
        for a, b in ((P_REG, P_INT), (P_INT, P_CORE)):
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=12,
                         color=via_col, lw=1.8, shrinkA=15, shrinkB=14, zorder=2))
        ax.add_patch(Circle(P_INT, 0.11, facecolor=C_INT, edgecolor="white", lw=1.1, zorder=4))
        ax.text(P_INT[0] + 0.17, P_INT[1], inter, ha="left", va="center", fontsize=8.0, color="#666", zorder=5)
    ax.add_patch(Circle(P_REG, 0.16, facecolor=C_REG, edgecolor="white", lw=1.3, zorder=4))
    ax.text(P_REG[0], P_REG[1] + 0.26, reg, ha="center", va="bottom", fontsize=9.5,
            fontweight="bold", color=C_REG, zorder=5)
    ax.add_patch(Circle(P_CORE, 0.18, facecolor=C_CORE, edgecolor="white", lw=1.3, zorder=4))
    ax.text(P_CORE[0], P_CORE[1] - 0.27, core, ha="center", va="top", fontsize=9.5,
            fontweight="bold", color=C_CORE, zorder=5)
    ax.set_xlim(-1.35, 1.6); ax.set_ylim(-1.35, 1.3); ax.set_aspect("equal"); ax.axis("off")


def _block(cases, mode, outname, ncols=4, show_legend=True):
    from matplotlib.lines import Line2D
    n = len(cases)
    nrows = int(np.ceil(n / ncols))
    extra = 0.7 if show_legend else 0.1
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.9, nrows * 3.0 + extra), squeeze=False)
    axes = axes.ravel()
    for ax, c in zip(axes, cases):
        _draw_case(ax, mode, c["reg"], c["core"], c.get("inter", ""))
    for ax in axes[n:]:
        ax.axis("off")
    if show_legend:
        leg = [Line2D([0], [0], color="#c1121f", lw=2.3, label="OmniPath relation"),
               Line2D([0], [0], color="#2c6fbb", lw=1.9, label="GRN relation"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor="#1f6fb2",
                      markeredgecolor="white", markersize=10, linestyle="None", label="regulator"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor="#bdbdbd",
                      markeredgecolor="white", markersize=8, linestyle="None", label="intermediate gene"),
               Line2D([0], [0], marker="o", color="none", markerfacecolor="#6a51a3",
                      markeredgecolor="white", markersize=11, linestyle="None", label="disease gene")]
        fig.legend(handles=leg, loc="lower center", ncol=5, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.01))
        fig.tight_layout(rect=[0, 0.06, 1, 1.0])
    else:
        fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/{outname}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/{outname}.pdf  ({n} cases)")


def fig_network(max_directed=6, max_via=6, max_ourvia=6):
    """Three exemplar-style figures (one per corroboration type), disease gene at the
    bottom of each subplot, with the intermediate shown when the relation is indirect."""
    e = pd.read_csv(f"{GR}/ALL.grn_edges.omnipath.tsv", sep="\t")
    pa = pd.read_csv(f"{GR}/ALL.paths_omnipath.tsv", sep="\t")

    # (1) directed in both: corroborated edge straight into a disease gene
    d = e[(e.op_directed) & (e.layer_dst == 0)][["reg_name", "tgt_name"]].drop_duplicates()
    cat1 = [{"reg": r.reg_name, "core": r.tgt_name} for r in d.head(max_directed).itertuples(index=False)]

    # (2) our direct, OmniPath via an intermediate (into a disease gene); show a few
    v = e[(e.op_via_middle) & (~e.op_directed) & (e.layer_dst == 0)][
        ["reg_name", "tgt_name", "via_middle_eg"]].drop_duplicates()
    cat2 = []
    for r in v.head(max_via).itertuples(index=False):
        eg = str(r.via_middle_eg).split(";")
        inter = eg[0] + (f" +{len(eg)-1}" if len(eg) > 1 and eg[0] else "")
        cat2.append({"reg": r.reg_name, "core": r.tgt_name, "inter": inter})

    # (3) our 2-hop via intermediate, OmniPath direct (the reference exemplar case)
    s = pa[pa.op_shortcut_GP_C][["grandparent", "parent", "child"]].drop_duplicates()
    s = s.groupby(["grandparent", "child"])["parent"].apply(list).reset_index().head(max_ourvia)
    cat3 = []
    for r in s.itertuples(index=False):
        inter = r.parent[0] + (f" +{len(r.parent)-1}" if len(r.parent) > 1 else "")
        cat3.append({"reg": r.grandparent, "core": r.child, "inter": inter})

    import matplotlib.gridspec as gridspec
    from matplotlib.lines import Line2D
    # transposed: one ROW per corroboration type, the cases laid out across columns
    rows = [("(a) direct in both", cat1, "both"),
            ("(b) GRN-direct,\nOmniPath via a gene", cat2, "opvia"),
            ("(c) GRN via a gene,\nOmniPath direct", cat3, "ourvia")]
    nrows = len(rows)
    ncases = max(len(c) for _, c, _ in rows)
    hr = [1.0] * nrows + [0.30]               # case rows + legend row
    wr = [0.95] + [1.0] * ncases              # row-label column + case columns
    fig = plt.figure(figsize=(sum(wr) * 2.3, sum(hr) * 2.6))
    gs = gridspec.GridSpec(len(hr), len(wr), figure=fig, height_ratios=hr,
                           width_ratios=wr, hspace=0.12, wspace=0.06)
    for i, (title, cases, mode) in enumerate(rows):
        lax = fig.add_subplot(gs[i, 0]); lax.axis("off")
        lax.text(0.96, 0.5, title, fontsize=12, fontweight="bold", ha="right", va="center")
        for j in range(ncases):
            ax = fig.add_subplot(gs[i, 1 + j])
            if j < len(cases):
                c = cases[j]
                _draw_case(ax, mode, c["reg"], c["core"], c.get("inter", ""))
            else:
                ax.axis("off")
    legax = fig.add_subplot(gs[-1, :]); legax.axis("off")
    leg = [Line2D([0], [0], color="#c1121f", lw=2.3, label="OmniPath relation"),
           Line2D([0], [0], color="#2c6fbb", lw=1.9, label="GRN relation"),
           Line2D([0], [0], marker="o", color="none", markerfacecolor="#1f6fb2",
                  markeredgecolor="white", markersize=10, linestyle="None", label="regulator"),
           Line2D([0], [0], marker="o", color="none", markerfacecolor="#bdbdbd",
                  markeredgecolor="white", markersize=8, linestyle="None", label="intermediate gene"),
           Line2D([0], [0], marker="o", color="none", markerfacecolor="#6a51a3",
                  markeredgecolor="white", markersize=11, linestyle="None", label="disease gene")]
    legax.legend(handles=leg, loc="center", ncol=5, fontsize=9.5, frameon=False)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/fig_flux_omnipath_network.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_flux_omnipath_network.pdf  "
          f"(a={len(cat1)} b={len(cat2)} c={len(cat3)})")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig_null()
    fig_network()
