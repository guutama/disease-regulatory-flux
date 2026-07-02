#!/usr/bin/env python3
"""Disease-regulatory flux-map figure of the manuscript (``fig_flux_map_grn``, ``fig:fluxgrn``).

The flux map drawn with the full GRN cascade: instead of drawing every regulator as a direct
edge into the disease gene, the two-hop structure is shown explicitly -- grandparent (hop 2) ->
parent (hop 1) -> disease gene. Edge width is |flux(ancestor -> gene)|, colour is the sign, and a
gold ring marks an ancestor that is itself a disease gene. The genes shown are the most
trans-driven, convergent (coherence > 0.7) gene-tissue cases pooled across the seven tissues.

Inputs (per-gene / per-edge SUMMARY tables only -- no individual-level data), under RESULTS:
  <AT>/association/association_<tissue>_<trait>.tsv   per-gene association (disease genes, z)
  <FL>/flux_edges_<tissue>_<trait>.tsv                regulator -> target flux edges
  <FEAT>/<tissue>.trans_features.tsv.gz               gene -> GRN-ancestor regulators (hop)
  <NET>/dag_<tissue>_orig_kde_fdr15_alltargets.csv    reconstructed GRN edges

Output: fig_flux_map_grn.{pdf,png} in the shared figures directory.

Run:  GENE_ANNOT=<gencode.v19.genes.tsv> python analysis/figures/plot_flux_map_grn.py

Input locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS); figures
are written to results/figures/ by default (override with FLUX_FIGURES); gene symbols come from
the annotation in $GENE_ANNOT. Style: Arial, 300 DPI.
"""
import sys, os, numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils")); import gene_labels as gl

TRAIT = "CAD_aragam2022"
# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS = Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_cv"))
FL = f"{RESULTS}/horseshoe_alltargets/flux"
AT = f"{RESULTS}/horseshoe_alltargets"
NET = f"{RESULTS}/network_alltargets"
FEAT = f"{RESULTS}/features_alltargets"
# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
OUT = str(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))
os.makedirs(OUT, exist_ok=True)
TISSUES = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
C_CORE, C_POS, C_NEG = "#6a51a3", "#c1121f", "#2c6fbb"
C_PARENT, C_GP, C_INT, C_DG = "#2171b5", "#9ecae1", "#cfcfcf", "#e8a000"
NTOP_GENES, NP, NG = 18, 5, 6
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "pdf.fonttype": 42, "savefig.dpi": 300})


def main():
    MAP = gl.load_label_map(gl.default_annot(os.environ.get("GENE_ANNOT")), "gene_id", "gene_name")
    data, cands = {}, []
    for t in TISSUES:
        A = pd.read_csv(f"{AT}/association/association_{t}_{TRAIT}.tsv", sep="\t")
        fe = pd.read_csv(f"{FL}/flux_edges_{t}_{TRAIT}.tsv", sep="\t")
        tf = pd.read_csv(f"{FEAT}/{t}.trans_features.tsv.gz", sep="\t",
                         usecols=["gene_id", "hop", "source_gene_id"]).drop_duplicates()
        hop = {(g, a): h for g, h, a in zip(tf.gene_id, tf.hop, tf.source_gene_id)}
        dag = pd.read_csv(f"{NET}/dag_{t}_orig_kde_fdr15_alltargets.csv", usecols=["Source", "Target"])
        dage = set(zip(dag.Source, dag.Target))
        par_of = dag.groupby("Target").Source.apply(set).to_dict()         # target -> {parents}
        data[t] = dict(fe=fe, hop=hop, dage=dage, par_of=par_of,
                       sig_genes=set(A.loc[A.p_adj < 0.05, "gene_id"]))
        grp = fe.groupby("target")
        nreg = grp.regulator.nunique(); conv = grp.flux.sum().abs() / grp.abs_flux.sum()
        s = A[(A.p_adj < 0.05) & A.z_trans.notna()].copy()
        s = s[(s.gene_id.map(nreg).fillna(0) >= 4) & (s.gene_id.map(conv).fillna(0) > 0.7)]
        for r in s.itertuples(index=False):
            cands.append({"tissue": t, "gene_id": r.gene_id, "z_trans": r.z_trans,
                          "z_twas": r.z_twas, "conv": float(conv.get(r.gene_id, np.nan))})
    C = pd.DataFrame(cands)
    C = C.reindex(C.z_trans.abs().sort_values(ascending=False).index).head(NTOP_GENES).reset_index(drop=True)

    # assemble per-gene cascades
    panels, gmax = [], 0.0
    for r in C.itertuples(index=False):
        t = r.tissue; fe = data[t]["fe"]; hop = data[t]["hop"]; dage = data[t]["dage"]
        sub = fe[fe.target == r.gene_id]
        pflux = {x.regulator: x.flux for x in sub.itertuples() if hop.get((r.gene_id, x.regulator)) == 1}
        gflux = {x.regulator: x.flux for x in sub.itertuples() if hop.get((r.gene_id, x.regulator)) == 2}
        parents = sorted(pflux, key=lambda p: -abs(pflux[p]))[:NP]
        gps = sorted(gflux, key=lambda g: -abs(gflux[g]))[:NG]
        par_all = data[t]["par_of"].get(r.gene_id, set())
        # link each grandparent to a routing parent (prefer a shown flux parent)
        link, inter = {}, {}
        for gp in gps:
            routing = [p for p in par_all if (gp, p) in dage]
            shown = [p for p in routing if p in parents]
            if shown:
                link[gp] = max(shown, key=lambda p: abs(pflux[p]))
            elif routing:
                link[gp] = routing[0]; inter[routing[0]] = True
            else:
                link[gp] = None
        panels.append((t, r.gene_id, r.z_trans, r.z_twas, r.conv, parents, gps,
                       pflux, gflux, link, list(inter)))
        gmax = max([gmax] + [abs(v) for v in pflux.values()] + [abs(v) for v in gflux.values()])

    ncols = 6; nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.0, nrows * 4.5)); axes = axes.ravel()
    YC, YP, YG = -1.45, -0.1, 1.25
    ew = lambda f: 0.6 + 4.6 * (abs(f) / gmax)

    for ax, (t, g, zt, ztot, cv, parents, gps, pflux, gflux, link, inter) in zip(axes, panels):
        sig_genes = data[t]["sig_genes"]
        midnodes = parents + [p for p in inter if p not in parents]
        px = {p: x for p, x in zip(midnodes, ([0.0] if len(midnodes) == 1
              else np.linspace(-1.15, 1.15, len(midnodes))))}
        core = np.array([0.0, YC])

        def node(p, xy, fill, sig_chk=True):
            if sig_chk and p in sig_genes:
                ax.add_patch(Circle(xy, 0.135, facecolor="none", edgecolor=C_DG, lw=2.2, zorder=2))
            ax.add_patch(Circle(xy, 0.08, facecolor=fill, edgecolor="white", lw=1.0, zorder=4))

        def edge(a, b, f, struct=False):
            col = "#bbbbbb" if struct else (C_POS if f >= 0 else C_NEG)
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=10,
                         lw=1.0 if struct else ew(f), color=col, alpha=0.6 if struct else 0.85,
                         shrinkA=9, shrinkB=12, zorder=2 if struct else 3))
            if not struct:
                d = b - a
                L = float(np.hypot(d[0], d[1])) or 1.0
                ang = np.degrees(np.arctan2(d[1], d[0]))
                if ang > 90:                                # keep the text upright/readable
                    ang -= 180
                elif ang < -90:
                    ang += 180
                perp = np.array([-d[1], d[0]]) / L          # unit vector across the arrow
                mp = a + 0.5 * d + 0.13 * perp              # sit beside the arrow, not on it
                ax.text(mp[0], mp[1], f"{f:+.2f}", ha="center", va="center", fontsize=7.4,
                        color=col, zorder=6, rotation=ang, rotation_mode="anchor",
                        bbox=dict(boxstyle="round,pad=0.06", fc="white", ec="none", alpha=0.7))

        # parents -> gene
        for p in midnodes:
            pp = np.array([px[p], YP])
            if p in pflux:
                edge(pp, core, pflux[p])
            else:
                edge(pp, core, 0, struct=True)                       # intermediate routing parent
            node(p, pp, C_PARENT if p in pflux else C_INT)
            ax.text(px[p], YP - 0.16, gl.resolve(str(p), MAP), ha="center", va="top",
                    fontsize=6.4, fontweight="bold" if p in sig_genes else "normal")
        # grandparents -> parent (or gene if unlinked); spread around their parent
        bypar = {}
        for gp in gps:
            bypar.setdefault(link[gp], []).append(gp)
        for tgt, glist in bypar.items():
            cx = px.get(tgt, 0.0)
            gxs = [cx] if len(glist) == 1 else np.linspace(cx - 0.5, cx + 0.5, len(glist))
            for gp, gx in zip(glist, gxs):
                gp_xy = np.array([gx, YG])
                dest = np.array([px[tgt], YP]) if tgt in px else core
                edge(gp_xy, dest, gflux[gp])
                node(gp, gp_xy, C_GP)
                ax.text(gx, YG + 0.12, gl.resolve(str(gp), MAP), ha="center", va="bottom",
                        fontsize=6.0, rotation=35, rotation_mode="anchor",
                        fontweight="bold" if gp in sig_genes else "normal")

        ax.add_patch(Circle(core, 0.155, facecolor=C_CORE, edgecolor="white", lw=1.4, zorder=5))
        ax.text(-0.04, YC - 0.24, f"{t} ", ha="right", va="top", fontsize=8.5, fontweight="bold", color="#888")
        ax.text(-0.02, YC - 0.24, gl.resolve(str(g), MAP), ha="left", va="top", fontsize=10.5,
                fontweight="bold", color=C_CORE)
        ax.text(0, 2.02, rf"$z_{{\mathrm{{total}}}}={ztot:+.1f}$   $z_{{\mathrm{{trans}}}}={zt:+.1f}$   "
                rf"conv $={cv:.2f}$", ha="center", va="top", fontsize=10.0, color="#333")
        ax.set_xlim(-1.75, 1.75); ax.set_ylim(-1.95, 2.05); ax.axis("off")
    for ax in axes[len(panels):]:
        ax.axis("off")

    leg = [Line2D([0], [0], marker="o", color="w", markerfacecolor=C_CORE, markersize=12, label="disease gene"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=C_PARENT, markersize=9, label="parent (hop 1)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=C_GP, markersize=8, label="grandparent (hop 2)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor=C_INT, markersize=8, label="intermediate (no flux)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor=C_DG,
                  markeredgewidth=2.2, markersize=13, label="ancestor also a disease gene"),
           Line2D([0], [0], color=C_POS, lw=3, label="flux toward disease risk"),
           Line2D([0], [0], color=C_NEG, lw=3, label="flux away from risk")]
    fig.legend(handles=leg, loc="lower center", ncol=7, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.07, hspace=0.18, wspace=0.06)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/fig_flux_map_grn.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_flux_map_grn.pdf  ({len(panels)} genes, "
          f"tissues {C.tissue.value_counts().to_dict()})")


if __name__ == "__main__":
    main()
