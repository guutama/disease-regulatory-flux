#!/usr/bin/env python3
"""Disease-regulatory flux-map cascade figure (fig_flux_map_grn): a grid of per-gene cascades
showing, for the most convergent disease genes, the hop-1 parents and hop-2 grandparents that
feed them and the flux each contributes.

Each panel: grandparent (hop 2) -> parent (hop 1) -> disease gene; edge width = |flux|, colour =
sign (red toward disease risk, blue away), a gold ring marks an ancestor that is itself a disease
gene. The displayed gene-tissue cases are the most convergent (coherence > 0.7 and >= 4
regulators), ranked by convergence x number of regulators and pooled across tissues.

Inputs (per tissue; tissues discovered from --flux-dir unless --tissues is given):
  --flux-dir    flux_<tissue>_<trait>.csv          S4 flux edges (regulator,...,target,...,hop,flux)
  --assoc-dir   association_<tissue>_<trait>.tsv    z_trans / z_twas / p_adj per gene
  --feat-dir    <tissue>.trans_features.tsv.gz      hop of each (gene, ancestor)
  --dag-dir     dag_<tissue>_<dag-tag>.csv          the GRN DAG (Source, Target) for parent links
  --gencode     GENCODE genes tsv (gene_id, gene_name) for HGNC symbols

Output (--outdir): fig_flux_map_grn.{pdf,png}

Usage:
  plot_flux_cascades.py --flux-dir DIR --assoc-dir DIR --feat-dir DIR --dag-dir DIR \
      --gencode FILE --outdir DIR [--trait CAD_aragam2022] [--dag-tag TAG] [--tissues ...]
"""
import sys, os, numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.lines import Line2D
import argparse, glob
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils")); import gene_labels as gl

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--flux-dir", required=True, help="dir with flux_<tissue>_<trait>.csv (S4 flux edges)")
ap.add_argument("--assoc-dir", required=True, help="dir with association_<tissue>_<trait>.tsv")
ap.add_argument("--feat-dir", required=True, help="dir with <tissue>.trans_features.tsv.gz")
ap.add_argument("--dag-dir", required=True, help="dir with the GRN DAG csv (dag_<tissue>_<dag-tag>.csv)")
ap.add_argument("--gencode", required=True, help="GENCODE genes tsv (gene_id, gene_name)")
ap.add_argument("--outdir", required=True, help="output dir for the figure")
ap.add_argument("--trait", default="CAD_aragam2022", help="trait label in the filenames")
ap.add_argument("--dag-tag", default="orig_kde_fdr15_alltargets", help="DAG filename tag: dag_<tissue>_<tag>.csv")
ap.add_argument("--tissues", nargs="+", default=None, help="tissue labels; default: discovered from --flux-dir")
args = ap.parse_args()

TRAIT = args.trait
FL = args.flux_dir              # flux_<t>_<trait>.csv (S4 columns)
AT = args.assoc_dir            # association_<t>_<trait>.tsv
NET = args.dag_dir             # dag_<t>_<dag-tag>.csv
FEAT = args.feat_dir           # <t>.trans_features.tsv.gz
ANNOT = args.gencode
OUT = args.outdir
DAG_TAG = args.dag_tag
os.makedirs(OUT, exist_ok=True)
TISSUES = args.tissues or sorted(
    os.path.basename(f)[len("flux_"):].rsplit(f"_{TRAIT}.csv", 1)[0]
    for f in glob.glob(f"{FL}/flux_*_{TRAIT}.csv"))
C_CORE, C_POS, C_NEG = "#6a51a3", "#c1121f", "#2c6fbb"
C_PARENT, C_GP, C_INT, C_DG = "#2171b5", "#9ecae1", "#cfcfcf", "#e8a000"
NCOLS = 6
NROWS = 6
NTOP_GENES, NP, NG = NCOLS * NROWS, 5, 6      # 6 rows x 6 cols = 36 gene-tissue cascade panels
CONV_MIN, NREG_MIN = 0.7, 4
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "pdf.fonttype": 42, "savefig.dpi": 300})


def read_flux(t):
    """S4-style flux CSV -> the figure's edge frame in Ensembl-id space with abs_flux.

    The figure operates on Ensembl ids (association/trans-features/DAG are all Ensembl), so the
    *_gene_id columns become `regulator`/`target`; the symbol columns are dropped."""
    fe = pd.read_csv(f"{FL}/flux_{t}_{TRAIT}.csv")
    # drop the HGNC-symbol columns; the figure works in Ensembl-id space
    fe = fe.drop(columns=["regulator", "target"])
    fe = fe.rename(columns={"regulator_gene_id": "regulator", "target_gene_id": "target"})
    fe["flux"] = fe["flux"].astype(float)
    fe["abs_flux"] = fe["flux"].abs()
    return fe[["regulator", "target", "flux", "abs_flux"]]


def main():
    MAP = gl.load_label_map(ANNOT, "gene_id", "gene_name")
    data, cands = {}, []
    for t in TISSUES:
        A = pd.read_csv(f"{AT}/association_{t}_{TRAIT}.tsv", sep="\t")
        fe = read_flux(t)
        tf = pd.read_csv(f"{FEAT}/{t}.trans_features.tsv.gz", sep="\t",
                         usecols=["gene_id", "hop", "source_gene_id"]).drop_duplicates()
        hop = {(g, a): h for g, h, a in zip(tf.gene_id, tf.hop, tf.source_gene_id)}
        dag = pd.read_csv(f"{NET}/dag_{t}_{DAG_TAG}.csv", usecols=["Source", "Target"])
        dage = set(zip(dag.Source, dag.Target))
        par_of = dag.groupby("Target").Source.apply(set).to_dict()         # target -> {parents}
        data[t] = dict(fe=fe, hop=hop, dage=dage, par_of=par_of,
                       sig_genes=set(A.loc[A.p_adj < 0.05, "gene_id"]))
        grp = fe.groupby("target")
        nreg = grp.regulator.nunique(); conv = grp.flux.sum().abs() / grp.abs_flux.sum()
        s = A[(A.p_adj < 0.05) & A.z_trans.notna()].copy()
        s = s[(s.gene_id.map(nreg).fillna(0) >= NREG_MIN) & (s.gene_id.map(conv).fillna(0) > CONV_MIN)]
        for r in s.itertuples(index=False):
            cands.append({"tissue": t, "gene_id": r.gene_id, "z_trans": r.z_trans,
                          "z_twas": r.z_twas, "conv": float(conv.get(r.gene_id, np.nan)),
                          "nreg": float(nreg.get(r.gene_id, 0))})
    C = pd.DataFrame(cands)
    print(f"[flux-map] candidate gene-tissue cases at conv>{CONV_MIN}, nreg>={NREG_MIN}: {len(C)}")
    if len(C) < NTOP_GENES:
        # relax the convergence floor just enough to fill the 6x6 grid
        need = NTOP_GENES
        relaxed = 0.5
        C2 = []
        for t in TISSUES:
            A = pd.read_csv(f"{AT}/association_{t}_{TRAIT}.tsv", sep="\t")
            fe = data[t]["fe"]; grp = fe.groupby("target")
            nreg = grp.regulator.nunique(); conv = grp.flux.sum().abs() / grp.abs_flux.sum()
            s = A[(A.p_adj < 0.05) & A.z_trans.notna()].copy()
            s = s[(s.gene_id.map(nreg).fillna(0) >= NREG_MIN) & (s.gene_id.map(conv).fillna(0) > relaxed)]
            for r in s.itertuples(index=False):
                C2.append({"tissue": t, "gene_id": r.gene_id, "z_trans": r.z_trans,
                           "z_twas": r.z_twas, "conv": float(conv.get(r.gene_id, np.nan)),
                           "nreg": float(nreg.get(r.gene_id, 0))})
        C = pd.DataFrame(C2)
        print(f"[flux-map] relaxed conv>{relaxed} to fill grid: {len(C)} cases")
    # rank by convergence AND regulator breadth (coherent convergence over many regulators),
    # not convergence alone
    C["score"] = C["conv"] * C["nreg"]
    C = C.sort_values("score", ascending=False).head(NTOP_GENES).reset_index(drop=True)

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

    ncols = NCOLS; nrows = int(np.ceil(len(panels) / ncols))
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
        # grandparents -> routing parent (or gene): each gene appears ONCE, spread to distinct
        # x across the top row (sorted by routing-parent x to keep edges tidy) so nodes/labels
        # never overlap. A gene already shown as a parent is not redrawn as a grandparent.
        gps_u = [gp for gp in gps if gp not in midnodes]
        gps_u = sorted(gps_u, key=lambda gp: (px.get(link[gp], 0.0), -abs(gflux[gp])))
        gxs = [0.0] if len(gps_u) <= 1 else np.linspace(-1.4, 1.4, len(gps_u))
        for gp, gx in zip(gps_u, gxs):
            gp_xy = np.array([gx, YG])
            tgt = link[gp]
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
    fig.legend(handles=leg, loc="lower center", ncol=7, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.01))
    fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.05, hspace=0.18, wspace=0.06)
    for ext in ("pdf", "png"):
        fig.savefig(f"{OUT}/fig_flux_map_grn.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}/fig_flux_map_grn.pdf  ({len(panels)} panels, {nrows} rows x {ncols} cols, "
          f"tissues {C.tissue.value_counts().to_dict()})")


if __name__ == "__main__":
    main()
