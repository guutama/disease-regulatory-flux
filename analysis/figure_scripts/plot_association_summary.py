#!/usr/bin/env python3
"""Combined figure for the disease-association step: one multi-panel summary of the
transcriptome-wide association of every gene's predicted expression, across tissues, plus a
stand-alone Manhattan.

  (a) |TWAS Z| by predictor class (cis-only / cis+trans / trans-only)
  (b) significant genes per tissue (FDR<0.05), stacked by class, trans-using % on top
  (c) |TWAS Z| vs the predictor SD sigma_g -- the total Z carries no sigma_g-driven inflation
      (reports r(|Z|,1/sigma_g)) and the residual trans-channel correlation among sig hits
  (d) transcriptome-wide Manhattan (best p per gene across tissues), rendered by the R
      ggplot2 + ggrepel helper (utils/manhattan_ggplot.R via plot_manhattan) when available,
      else a matplotlib fallback

Inputs (per tissue; tissues are discovered from --assoc-dir unless --tissues is given):
  --assoc-dir   association_<tissue>_<trait>.tsv   (from the TWAS association step)
  --gencode     GENCODE genes tsv, for the Manhattan genomic positions

Output (--outdir):
  fig_assoc_combined.{pdf,png}               the multi-panel figure
  fig_assoc_manhattan_alltargets.{pdf,png}   the stand-alone Manhattan

Usage:
  plot_association_summary.py --assoc-dir DIR --outdir DIR --gencode FILE \
      [--trait CAD_aragam2022] [--tissues T1 T2 ...]
"""
import argparse
import glob
import os
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import norm

# shared plotting helpers live in <repo>/utils; resolve relative to this file (CWD-independent)
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))
import gene_labels as gl                      # noqa: F401  (shared label helpers)
import plot_manhattan as M

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--assoc-dir", required=True,
                help="dir with association_<tissue>_<trait>.tsv from the association step")
ap.add_argument("--outdir", required=True, help="output dir for the figure(s)")
ap.add_argument("--gencode", required=True,
                help="GENCODE genes tsv (chr, gene_id, gene_name, start, end, ...) for Manhattan positions")
ap.add_argument("--trait", default="CAD_aragam2022",
                help="trait label in the association filenames (default: CAD_aragam2022)")
ap.add_argument("--tissues", nargs="+", default=None,
                help="tissue labels in plotting order; default: discovered from --assoc-dir")
args = ap.parse_args()

TRAIT = args.trait
AD = args.assoc_dir
OUT = args.outdir
ANNOT = args.gencode
TIS = args.tissues or sorted(
    os.path.basename(f)[len("association_"):].rsplit(f"_{TRAIT}.tsv", 1)[0]
    for f in glob.glob(f"{AD}/association_*_{TRAIT}.tsv"))
if not TIS:
    raise SystemExit(f"no association_*_{TRAIT}.tsv found in {AD}")
C_CIS, C_BOTH, C_TRANS = "#4C72B0", "#55A868", "#DD8452"
LAB = {"cis_only": "cis only", "cis_trans": "cis+trans", "trans_only": "trans only"}
COL = {"cis_only": C_CIS, "cis_trans": C_BOTH, "trans_only": C_TRANS}
NULL = (2 / np.pi) ** 0.5
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 9, "axes.linewidth": 0.8, "savefig.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False})

os.makedirs(OUT, exist_ok=True)
A = pd.concat([pd.read_csv(f"{AD}/association_{t}_{TRAIT}.tsv", sep="\t") for t in TIS], ignore_index=True)
# Rename columns to what the plotting code below expects: p_value -> p, config -> win_config.
A["p"] = A["p_value"]
A["win_config"] = A["config"]
sig = A[A["p_adj"] < 0.05].copy()
zcut = abs(norm.isf(sig["p"].max() / 2))

# ---- (d) Manhattan: use the R ggrepel renderer if available, else the matplotlib fallback ----
M.ASSOC_DIR = Path(AD)
M.ASSOC_FMT = "association_{tissue}_" + TRAIT + ".tsv"
M.GENCODE = Path(ANNOT)
M.TITLE = "Transcriptome-wide association with coronary artery disease"
# best (min-p) row per unique gene; build it here and reuse M only for genomic positions.
best = A.sort_values("p").drop_duplicates(subset="gene_id").reset_index(drop=True)
best = M.attach_positions(best, M.load_gene_positions())
n_pos = int(best["chrom"].notna().sum())
n_sig_mapped = int((best["p_adj"] < M.FDR).sum())
print(f"manhattan: {len(best)} unique genes, {n_pos} with GENCODE position "
      f"({100*n_pos/len(best):.1f}% coverage), {n_sig_mapped} FDR<{M.FDR}")
mh_pdf = Path(OUT) / "fig_assoc_manhattan_alltargets.pdf"
ok = M.try_ggplot_manhattan(best, mh_pdf, M.LABEL_MIN_NEGLOG, 50.0, width=18, height=6)
mh_png = mh_pdf.with_suffix(".png")

fig = plt.figure(figsize=(17, 11))
sf_top, sf_bot = fig.subfigures(2, 1, height_ratios=[0.8, 1.0])
axA, axB, axC = sf_top.subplots(1, 3, gridspec_kw={"width_ratios": [1, 2.0, 1.25]})
sf_top.subplots_adjust(left=0.045, right=0.99, top=0.84, bottom=0.20, wspace=0.22)

# ---- (a) |Z| by class ----
order = ["cis_only", "cis_trans", "trans_only"]
data = [A.loc[A["win_config"] == c, "z_twas"].abs().dropna().values for c in order]
vp = axA.violinplot(data, showextrema=False, widths=0.85)
for b, c in zip(vp["bodies"], order):
    b.set_facecolor(COL[c]); b.set_alpha(0.75); b.set_edgecolor("white")
for i, d in enumerate(data, 1):
    axA.hlines(np.median(d), i - 0.25, i + 0.25, color="#222", lw=1.4, zorder=3)
axA.axhline(NULL, ls="--", lw=0.9, color="#888")
axA.axhline(zcut, ls=":", lw=0.9, color="#c0392b")
axA.set_xticks([1, 2, 3]); axA.set_xticklabels([LAB[c] for c in order], rotation=12)
axA.set_ylabel("|TWAS $Z$|"); axA.set_ylim(0, np.percentile(np.concatenate(data), 99.5))
axA.set_title("(a) TWAS signal by class", loc="left", fontweight="bold")
axA.text(0.97, 0.97, "··· FDR 5%\n– – null $E|Z|$", transform=axA.transAxes,
         ha="right", va="top", fontsize=6.8, color="#555")

# ---- (b) sig genes per tissue, stacked by class ----
ct = (sig.groupby(["tissue", "win_config"]).size().unstack(fill_value=0).reindex(TIS)[order])
x = np.arange(len(TIS)); bottom = np.zeros(len(TIS))
for c in order:
    axB.bar(x, ct[c].values, bottom=bottom, color=COL[c], width=0.72, label=LAB[c])
    bottom += ct[c].values
for i in range(len(TIS)):
    tot = bottom[i]; tu = ct.iloc[i][["cis_trans", "trans_only"]].sum()
    axB.text(i, tot + 8, f"{100 * tu / tot:.0f}%", ha="center", va="bottom", fontsize=7.5, color="#333")
axB.set_xticks(x); axB.set_xticklabels(TIS)
axB.set_ylabel("significant genes (FDR$<$0.05)")
axB.set_title("(b) Significant genes per tissue", loc="left", fontweight="bold")
axB.legend(frameon=False, ncol=3, fontsize=8, loc="upper right")

# ---- (c) |Z| vs predictor SD: large |Z| is not bought by tiny sigma_g ----
az = A["z_twas"].abs().values
axC.scatter(A["sigma_g"], az, s=4, c="#d9d9d9", alpha=0.45, linewidths=0,
            rasterized=True)                                   # all predictable genes
for c in order:                                                # significant, by class
    m = (A["win_config"] == c) & (A["p_adj"] < 0.05)
    axC.scatter(A.loc[m, "sigma_g"], A.loc[m, "z_twas"].abs(), s=7, c=COL[c],
                alpha=0.8, linewidths=0)
axC.axhline(zcut, ls=":", lw=0.9, color="#c0392b")
r = np.corrcoef(az, 1.0 / A["sigma_g"].values)[0, 1]
st = A[(A["p_adj"] < 0.05) & (A["n_trans"] > 0) & (A["sigma_g"] > 0.2)]   # significant trans-using genes above the trans sigma_g floor
rt = np.corrcoef(st["z_trans"].abs(), 1.0 / st["sigma_g"])[0, 1]
axC.text(0.96, 0.96,
         f"$r(|Z_{{\\mathrm{{twas}}}}|,\\,1/\\sigma_g)={r:.2f}$\n"
         f"$r(|Z_{{\\mathrm{{trans}}}}|,\\,1/\\sigma_g)={rt:.2f}$  (sig, $\\sigma_g>0.2$)",
         transform=axC.transAxes, ha="right", va="top", fontsize=7.5, color="#333")
axC.set_xlabel("$\\sigma_g$  (predicted-expression SD)")
axC.set_ylabel("$|Z_{\\mathrm{twas}}|$")
axC.set_title("(c) No $\\sigma_g$-driven inflation", loc="left", fontweight="bold")

# ---- (d) Manhattan full-width: R render if available, else matplotlib fallback ----
if ok and mh_png.exists():
    axc = sf_bot.subplots(1, 1)
    sf_bot.subplots_adjust(left=0.05, right=0.99, top=0.99, bottom=0.01)
    axc.imshow(mpimg.imread(mh_png)); axc.axis("off")
    axc.text(0.0, 1.0, "(d)", transform=axc.transAxes, fontsize=13, fontweight="bold",
             va="bottom", ha="left")
else:
    # Rscript / ggplot2 unavailable -> native matplotlib Manhattan (plot_manhattan.render):
    # best p per gene across the 7 tissues, genomic x from GENCODE, FDR threshold line, labels.
    sf_bot.subplots_adjust(left=0.06, right=0.99, top=0.90, bottom=0.11, hspace=0.06)
    M.render(sf_bot, best)
    sf_bot.text(0.008, 0.99, "(d)", fontsize=13, fontweight="bold", va="top", ha="left")

for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/fig_assoc_combined.{ext}", bbox_inches="tight")
print("wrote", f"{OUT}/fig_assoc_combined.pdf  (manhattan R render ok={ok})")
print(f"sig gene-tissue={len(sig)} unique={sig.gene_id.nunique()} "
      f"trans-only(unique)={sig[sig.win_config=='trans_only'].gene_id.nunique()}")
