#!/usr/bin/env python3
"""Background co-expression versus reconstructed-network density per tissue
(``fig_coexpression_density``).

Panel (a): the mean absolute pairwise expression correlation per tissue, computed two ways ---
on each tissue's own expressed genes and on the set of genes common to all tissues --- so the
comparison is not confounded by different gene sets. Both use equal-size random gene subsamples
(default 1,500 genes, averaged over several draws), correlated across the tissue's own samples.

Panel (b): network density (mean regulator out-degree, from the reconstructed DAG) against the
common-gene background co-expression, with the across-tissue Pearson correlation.

The two adipose depots (SF, VAF) are highlighted.

Inputs: the reconstructed per-tissue DAG edge lists under RESULTS (network density) and the
per-tissue expression matrices under --expr-dir (only aggregate co-expression is reported, never
individual values). Output: fig_coexpression_density.{pdf,png} in the shared figures directory.

Network inputs resolve under RESULTS (default results_findr_py/, override with FLUX_RESULTS); figures go
to results/figures/ by default (override with FLUX_FIGURES). Style: Arial, 300 DPI.

Usage:
    plot_coexpression_density.py [--net-dir DIR] [--expr-dir DIR] [--n-genes 1500]
                                 [--repeats 10] [--seed 0] [--out FILE] [TISSUE ...]
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8, "savefig.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False,
})

C_OWN, C_COMMON, C_ADIPOSE, C_OTHER = "#B0B7C0", "#4C72B0", "#C0392B", "#4C72B0"
ADIPOSE = {"SF", "VAF"}
# Input result locations resolve under RESULTS (default results_findr_py/, override with FLUX_RESULTS).
RESULTS = Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_findr_py"))
# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
FIGDIR = Path(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))


def load_expression(path: str):
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    with op as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi = header.index("gene_id")
        ncol = len(header) - 2
        genes, rows = [], []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genes.append(f[gi])
            rows.append(f[2:2 + ncol])
    return genes, np.array(rows, dtype=np.float64)


def mean_abs_corr(mat, gene_idx, n_genes, repeats, rng):
    vals = []
    for _ in range(repeats):
        pick = rng.choice(gene_idx, size=min(n_genes, gene_idx.size), replace=False)
        c = np.corrcoef(mat[pick])
        iu = np.triu_indices(c.shape[0], k=1)
        vals.append(np.nanmean(np.abs(c[iu])))
    return float(np.mean(vals))


def edges_and_outdegree(dag_path):
    sources, n = [], 0
    with open(dag_path) as fh:
        head = fh.readline().rstrip("\n").split(",")
        si = head.index("Source")
        for line in fh:
            sources.append(line.split(",", si + 1)[si]); n += 1
    n_reg = len(set(sources))
    return n, (n / n_reg if n_reg else 0.0)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("tissues", nargs="*",
                   default=["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"])
    p.add_argument("--net-dir", default=str(RESULTS / "network_alltargets"))
    p.add_argument("--expr-dir", default="results/stage2_test/expression")
    p.add_argument("--n-genes", type=int, default=1500)
    p.add_argument("--repeats", type=int, default=10)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default=str(FIGDIR / "fig_coexpression_density.pdf"))
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    expr, common = {}, None
    for t in args.tissues:
        genes, mat = load_expression(os.path.join(args.expr_dir, f"{t}.expr.tsv.gz"))
        expr[t] = (genes, mat)
        common = set(genes) if common is None else (common & set(genes))
    common_sorted = sorted(common)
    sys.stderr.write(f"[fig] {len(common_sorted)} common genes across {len(args.tissues)} tissues\n")

    rng = np.random.default_rng(args.seed)
    rec = []
    for t in args.tissues:
        genes, mat = expr[t]
        pos = {g: i for i, g in enumerate(genes)}
        own_idx = np.arange(len(genes))
        com_idx = np.array([pos[g] for g in common_sorted])
        mac_own = mean_abs_corr(mat, own_idx, args.n_genes, args.repeats, rng)
        mac_com = mean_abs_corr(mat, com_idx, args.n_genes, args.repeats, rng)
        n_edges, outdeg = edges_and_outdegree(
            os.path.join(args.net_dir, f"dag_{t}_orig_kde_fdr15_alltargets.csv"))
        rec.append(dict(tissue=t, edges=n_edges, outdeg=outdeg,
                        mac_own=mac_own, mac_com=mac_com))
        sys.stderr.write(f"  {t}: edges={n_edges} outdeg={outdeg:.1f} "
                         f"|corr|_own={mac_own:.3f} |corr|_common={mac_com:.3f}\n")

    rec.sort(key=lambda r: r["outdeg"])              # order tissues by network density
    tis = [r["tissue"] for r in rec]
    x = np.arange(len(tis))

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(8.2, 3.4))

    # (a) grouped bars: own-gene vs common-gene background co-expression
    w = 0.38
    axa.bar(x - w / 2, [r["mac_own"] for r in rec], w, color=C_OWN,
            edgecolor="black", linewidth=0.4, label="all expressed genes")
    axa.bar(x + w / 2, [r["mac_com"] for r in rec], w, color=C_COMMON,
            edgecolor="black", linewidth=0.4, label="common gene set")
    axa.set_xticks(x); axa.set_xticklabels(tis, rotation=0)
    axa.set_ylabel("mean |pairwise expression correlation|")
    axa.set_title("(a) Background co-expression per tissue", loc="left", fontsize=9.5)
    axa.legend(frameon=False, fontsize=8, loc="upper left")
    for i, r in enumerate(rec):
        if r["tissue"] in ADIPOSE:
            axa.annotate("adipose", (i, max(r["mac_own"], r["mac_com"])),
                         textcoords="offset points", xytext=(0, 3), ha="center",
                         fontsize=7, color=C_ADIPOSE, fontweight="bold")

    # (b) density vs common-gene background co-expression
    od = np.array([r["outdeg"] for r in rec])
    mac = np.array([r["mac_com"] for r in rec])
    colors = [C_ADIPOSE if r["tissue"] in ADIPOSE else C_OTHER for r in rec]
    axb.scatter(od, mac, c=colors, s=42, edgecolor="black", linewidth=0.5, zorder=3)
    b, a = np.polyfit(od, mac, 1)
    xs = np.linspace(od.min(), od.max(), 50)
    axb.plot(xs, a + b * xs, "--", color="grey", linewidth=1, zorder=1)
    label_off = {"LIV": (5, 4), "MAM": (5, -11), "Blood": (6, 4),
                 "AOR": (5, 4), "SKLM": (5, 4), "SF": (-6, 6), "VAF": (-8, 6)}
    for r in rec:
        dx, dy = label_off.get(r["tissue"], (4, 3))
        ha = "right" if dx < 0 else "left"
        axb.annotate(r["tissue"], (r["outdeg"], r["mac_com"]),
                     textcoords="offset points", xytext=(dx, dy), fontsize=7.5, ha=ha,
                     color=(C_ADIPOSE if r["tissue"] in ADIPOSE else "black"))
    axb.set_xlabel("mean regulator out-degree (network density)")
    axb.set_ylabel("mean |corr| (common gene set)")
    axb.set_title("(b) Co-expression tracks density", loc="left", fontsize=9.5)
    axb.text(0.97, 0.05, f"Pearson $r$ = {np.corrcoef(od, mac)[0, 1]:.2f}",
             transform=axb.transAxes, ha="right", va="bottom", fontsize=8.5)

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, bbox_inches="tight")
    fig.savefig(args.out.replace(".pdf", ".png"), bbox_inches="tight", dpi=300)
    sys.stderr.write(f"[fig] wrote {args.out}\n")


if __name__ == "__main__":
    main()
