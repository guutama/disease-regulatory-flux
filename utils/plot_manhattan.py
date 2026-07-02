#!/usr/bin/env python3
"""
Transcriptome-wide TWAS Manhattan plot for CAD (fig8a_manhattan_combined).

Best p-value per gene across the 7 STARNET tissues. Genes with -log10(p) above
--high-thresh (default 15) are labelled by HGNC symbol.

Primary renderer = ggplot2 + ggrepel via utils/manhattan_ggplot.R (the same
function used in EXP-GRN-TWAS: ggrepel spaces the labels so they do not overlap;
-log10(p) is capped at --cap, default 50, so one extreme peak does not compress
the rest). Falls back to a self-contained matplotlib renderer (with an automatic
broken y-axis) only if the R toolchain is unavailable.

Usage:
    python3 utils/plot_manhattan.py [--high-thresh 15] [--cap 50] [--no-r]
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_BIN = Path(__file__).resolve().parent
_ROOT = _BIN.parent
import sys
sys.path.insert(0, str(_ROOT / "utils")); sys.path.insert(0, str(_BIN))
try:
    import figure_data as fd
    fd.apply_style(base=8)
except Exception as e:
    print(f"  [style] figure_data.apply_style unavailable ({e}); using fallback rcParams")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    })

TISSUES = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
AUTOSOMES = list(range(1, 23))
FDR = 0.05

ASSOC_DIR = _ROOT / "results" / "stage10_nf" / "association"
ASSOC_FMT = "association_{tissue}_CAD_aragam2022.tsv"
GENCODE = Path(os.environ.get("GENE_ANNOT", ""))   # gene annotation; set via $GENE_ANNOT or overridden by the caller
OUT_DIR = _ROOT / "results" / "figures"
OUT_NAME = "fig8a_manhattan_combined"

RSCRIPT_WRAPPER = _BIN / "rscript_wrapper.sh"
MANHATTAN_R = _BIN / "manhattan_ggplot.R"
TITLE = "Transcriptome-wide association with coronary artery disease"

CHR_COLS = ["#4E79A7", "#A0CBE8"]
RED = "#CC0000"


def load_gene_positions():
    """gene_id (versioned ENSG) -> chrom (int), pos (midpoint). Also a base-ENSG map."""
    gc = pd.read_csv(GENCODE, sep="\t",
                     usecols=["chr", "gene_id", "start", "end", "gene_name"],
                     dtype={"chr": str})
    gc = gc[gc["chr"].isin([str(c) for c in AUTOSOMES])].copy()
    gc["chrom"] = gc["chr"].astype(int)
    gc["pos"] = ((gc["start"] + gc["end"]) // 2).astype(int)
    gc["ensg_base"] = gc["gene_id"].str.split(".").str[0]
    return gc[["gene_id", "ensg_base", "chrom", "pos", "gene_name"]]


def load_best_per_gene():
    frames = []
    for t in TISSUES:
        p = ASSOC_DIR / ASSOC_FMT.format(tissue=t)
        if not p.exists():
            print(f"  [skip] {p.name} not found")
            continue
        head = pd.read_csv(p, sep="\t", nrows=0).columns
        cols = [c for c in ("gene", "gene_id", "p", "p_adj") if c in head]
        df = pd.read_csv(p, sep="\t", usecols=cols)
        df["tissue"] = t
        frames.append(df)
    if not frames:
        raise SystemExit("No association files found.")
    assoc = pd.concat(frames, ignore_index=True)
    print(f"  {len(assoc):,} gene-tissue rows across {len(frames)} tissues")
    best = (assoc.sort_values("p")
                 .drop_duplicates(subset="gene_id")
                 .reset_index(drop=True))
    return best


def attach_positions(best, gc):
    # Direct match on versioned gene_id.
    pos = gc[["gene_id", "chrom", "pos", "gene_name"]]
    merged = best.merge(pos, on="gene_id", how="left")
    # If the association table had no HGNC 'gene' column, use GENCODE gene_name;
    # otherwise keep 'gene' but backfill any blanks from gene_name.
    if "gene" not in merged.columns:
        merged["gene"] = merged["gene_name"]
    else:
        merged["gene"] = merged["gene"].fillna(merged["gene_name"])
    merged["gene"] = merged["gene"].fillna(merged["gene_id"].str.split(".").str[0])
    # Fallback: version-stripped base ENSG for any misses.
    miss = merged["chrom"].isna()
    if miss.any():
        base_map = (gc.sort_values("gene_id")
                      .drop_duplicates("ensg_base")
                      .set_index("ensg_base")[["chrom", "pos"]])
        merged["ensg_base"] = merged["gene_id"].str.split(".").str[0]
        fill = merged.loc[miss, "ensg_base"].map(
            lambda b: base_map.loc[b].to_dict() if b in base_map.index else None)
        for idx, val in fill.items():
            if val is not None:
                merged.at[idx, "chrom"] = val["chrom"]
                merged.at[idx, "pos"] = val["pos"]
        merged = merged.drop(columns=["ensg_base"])
    return merged


# A broken y-axis is used only when one (or few) peaks sit far above the rest,
# so the bulk of the genes are not visually compressed. The split is triggered
# only when the gap between the top point and the next is large (data-agnostic).
BREAK_MIN_GAP = 40.0       # -log10(p) units between top cluster and the rest
BREAK_PAD = 8.0            # half-height of the upper (peak) panel
LABEL_MIN_NEGLOG = 10.0    # only label genes with -log10(p) above this


def _compute_xy(df):
    """Add neglog + cumulative genomic x_plot; return (df, chrom_mids, gap)."""
    df = df.dropna(subset=["chrom", "pos", "p"]).copy()
    df["chrom"] = df["chrom"].astype(int)
    df["pos"] = df["pos"].astype(int)
    df["neglog"] = -np.log10(df["p"].clip(lower=1e-300))
    gap = int(20e6)
    offsets, cum = {}, 0
    for c in AUTOSOMES:
        offsets[c] = cum
        sub_pos = df.loc[df["chrom"] == c, "pos"]
        cum += (int(sub_pos.max()) if len(sub_pos) else 0) + gap
    df["x_plot"] = df["chrom"].map(offsets).astype(float) + df["pos"]
    chrom_mids = {}
    for c in AUTOSOMES:
        sub = df.loc[df["chrom"] == c, "x_plot"]
        if len(sub):
            chrom_mids[c] = (sub.min() + sub.max()) / 2
    return df, chrom_mids, gap


def _scatter(ax, df):
    """Per-chromosome scatter (non-sig faint/small, sig red/larger on top)."""
    for c in AUTOSOMES:
        sub = df[df["chrom"] == c]
        if sub.empty:
            continue
        cidx = c % 2
        ns = sub[sub["p_adj"] >= FDR]
        sg = sub[sub["p_adj"] < FDR]
        ax.scatter(ns["x_plot"], ns["neglog"],
                   s=1.5, color=CHR_COLS[cidx], alpha=0.45, rasterized=True)
        ax.scatter(sg["x_plot"], sg["neglog"],
                   s=6, color=RED, alpha=0.9, zorder=3, rasterized=True)


def _label_sig(ax, df, only_above=None, only_below=None):
    """Annotate FDR-significant genes by HGNC symbol; clip_on keeps off-range
    labels out of each broken-axis panel. Optional neglog window filters which
    genes are labelled on a given panel. Only genes with -log10(p) above
    LABEL_MIN_NEGLOG are labelled (keeps the dense baseline uncluttered)."""
    sig = df[(df["p_adj"] < FDR) & (df["neglog"] > LABEL_MIN_NEGLOG)]
    if only_above is not None:
        sig = sig[sig["neglog"] >= only_above]
    if only_below is not None:
        sig = sig[sig["neglog"] < only_below]
    for _, row in sig.iterrows():
        ax.annotate(
            str(row["gene"]),
            xy=(row["x_plot"], row["neglog"]),
            xytext=(0, 4), textcoords="offset points",
            fontsize=5, ha="center", va="bottom", clip_on=True,
            arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4,
                            shrinkA=0, shrinkB=1),
        )


def render(fig, best):
    """Manhattan with an automatic broken y-axis for an isolated top peak."""
    df, chrom_mids, gap = _compute_xy(best)
    sig_mask = df["p_adj"] < FDR
    thresh_y = df.loc[sig_mask, "neglog"].min() if sig_mask.any() else -np.log10(FDR)
    xmin, xmax = df["x_plot"].min() - gap, df["x_plot"].max() + gap
    xticks = [chrom_mids[c] for c in sorted(chrom_mids)]
    xlabels = [str(c) if c % 2 == 1 else "" for c in sorted(chrom_mids)]

    order = np.sort(df["neglog"].values)[::-1]
    top_val = float(order[0])
    second = float(order[1]) if len(order) > 1 else top_val
    broken = (top_val - second) > BREAK_MIN_GAP

    title = (f"Transcriptome-wide association with coronary artery disease  "
             f"({int(sig_mask.sum()):,} significant genes)")

    if not broken:
        ax = fig.subplots(1, 1)
        _scatter(ax, df)
        ax.axhline(thresh_y, color=RED, lw=0.7, ls="--", alpha=0.7)
        _label_sig(ax, df)
        ax.set_xticks(xticks); ax.set_xticklabels(xlabels, fontsize=7)
        ax.set_xlim(xmin, xmax); ax.margins(y=0.08)
        ax.set_xlabel("Chromosome"); ax.set_ylabel(r"$-\log_{10}(p)$")
        ax.set_title(title, fontweight="bold")
        return

    # Broken axis: small upper panel for the isolated peak(s), large lower panel.
    ax_top, ax_bot = fig.subplots(
        2, 1, sharex=True, gridspec_kw=dict(height_ratios=[1, 9], hspace=0.06))
    cut = second + 0.10 * (top_val - second)        # boundary between panels
    _scatter(ax_top, df); _scatter(ax_bot, df)

    ax_top.set_ylim(top_val - BREAK_PAD, top_val + BREAK_PAD)
    ax_bot.set_ylim(0, np.ceil(cut / 5) * 5)
    ax_bot.axhline(thresh_y, color=RED, lw=0.7, ls="--", alpha=0.7)
    _label_sig(ax_bot, df, only_below=cut)
    _label_sig(ax_top, df, only_above=cut)

    # Hide the inner spines and draw the diagonal break marks.
    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_top.set_yticks([round(top_val)])
    d = 0.008
    for ax, ys in ((ax_top, (0, 0)), (ax_bot, (1, 1))):
        kw = dict(transform=ax.transAxes, color="k", clip_on=False, lw=0.8)
        sign = -1 if ax is ax_top else +1
        y = ys[0]
        ax.plot((-d, +d), (y + sign * d * 9, y - sign * d * 9), **kw)

    ax_bot.set_xticks(xticks); ax_bot.set_xticklabels(xlabels, fontsize=7)
    ax_bot.set_xlim(xmin, xmax)
    ax_bot.set_xlabel("Chromosome")
    fig.supylabel(r"$-\log_{10}(p)$", x=0.065, fontsize=10)
    ax_top.set_title(title, fontweight="bold")


def try_ggplot_manhattan(best, pdf_path, high_thresh, cap, width=18, height=6):
    """Render via ggplot2 + ggrepel (utils/manhattan_ggplot.R). Returns True on
    success. The R script writes both <out>.pdf and the sibling <out>.png."""
    if not (RSCRIPT_WRAPPER.exists() and MANHATTAN_R.exists()):
        print("  [ggplot2] R wrapper/script not found")
        return False
    plot_df = best.dropna(subset=["chrom", "pos", "p"]).copy()
    plot_df["chrom"] = plot_df["chrom"].astype(int)
    plot_df["pos"] = plot_df["pos"].astype(int)
    plot_df = plot_df.rename(columns={"gene": "gene_name", "p": "p_twas",
                                      "p_adj": "p_twas_adj"})
    cols = ["gene_name", "chrom", "pos", "p_twas", "p_twas_adj"]
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as tmp:
        tmp_csv = tmp.name
    try:
        plot_df[cols].to_csv(tmp_csv, index=False)
        cmd = ["bash", str(RSCRIPT_WRAPPER), str(MANHATTAN_R),
               "--input", tmp_csv, "--output", str(pdf_path),
               "--title", TITLE, "--fdr", str(FDR), "--color", "#CC0000",
               "--cap", str(cap), "--high-thresh", str(high_thresh),
               "--label-n", "0", "--width", str(width), "--height", str(height)]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        if res.returncode == 0 and pdf_path.exists():
            print("  " + (res.stdout.strip() or f"Saved {pdf_path}"))
            return True
        print(f"  [ggplot2] R failed (rc={res.returncode}):\n{res.stderr[-1500:]}")
        return False
    except Exception as e:
        print(f"  [ggplot2] subprocess error: {e}")
        return False
    finally:
        Path(tmp_csv).unlink(missing_ok=True)


def main():
    global LABEL_MIN_NEGLOG, ASSOC_DIR, OUT_DIR, OUT_NAME
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--high-thresh", type=float, default=LABEL_MIN_NEGLOG,
                    help="label genes with -log10(p) above this (default 15)")
    ap.add_argument("--cap", type=float, default=50.0,
                    help="-log10(p) ceiling for the ggrepel renderer (default 50)")
    ap.add_argument("--no-r", action="store_true",
                    help="skip the ggplot2/ggrepel renderer; use matplotlib")
    ap.add_argument("--assoc-dir", default=str(ASSOC_DIR),
                    help="dir of association_<tissue>_CAD_aragam2022.tsv tables")
    ap.add_argument("--out-dir", default=str(OUT_DIR), help="output directory")
    ap.add_argument("--out-name", default=OUT_NAME, help="output file basename")
    args = ap.parse_args()
    ASSOC_DIR = Path(args.assoc_dir); OUT_DIR = Path(args.out_dir); OUT_NAME = args.out_name

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LABEL_MIN_NEGLOG = args.high_thresh   # used by the matplotlib fallback labels

    print("Loading gene positions (GENCODE v19)...")
    gc = load_gene_positions()
    print(f"  {len(gc):,} autosomal gencode genes")

    print("Loading association results...")
    best = load_best_per_gene()
    print(f"  {len(best):,} unique genes (best p per gene_id)")

    best = attach_positions(best, gc)
    print(f"  {int(best['chrom'].notna().sum()):,} genes with genomic position")
    print(f"  FDR < {FDR}: {int((best['p_adj'] < FDR).sum()):,} significant genes")
    print(f"  labelled (-log10 p > {args.high_thresh:g}): "
          f"{int((-np.log10(best['p'].clip(lower=1e-300)) > args.high_thresh).sum()):,}")

    pdf_path = OUT_DIR / f"{OUT_NAME}.pdf"
    png_path = OUT_DIR / f"{OUT_NAME}.png"

    if not args.no_r:
        print("Rendering with ggplot2 + ggrepel...")
        if try_ggplot_manhattan(best, pdf_path, args.high_thresh, args.cap):
            return
        print("  -> falling back to matplotlib renderer")

    fig = plt.figure(figsize=(14, 5))
    render(fig, best)
    fig.savefig(pdf_path, bbox_inches="tight", dpi=200)
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  Saved {pdf_path}")
    print(f"  Saved {png_path}")


if __name__ == "__main__":
    main()
