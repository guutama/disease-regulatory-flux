#!/usr/bin/env python3
"""Combined CAD transcriptome-wide association figure of the manuscript (``fig_assoc_combined``).

Three panels consolidating the association results:
  (a) |TWAS Z| by predictor class (cis / cis+trans / trans-only);
  (b) significant genes per tissue, stacked by class (trans-using % on top);
  (c) a transcriptome-wide Manhattan (best p per gene across the seven tissues), rendered by the
      ggplot2 + ggrepel script utils/manhattan_ggplot.R via plot_manhattan.try_ggplot_manhattan
      and composited full-width. If R is unavailable, a matplotlib Manhattan is used instead.

Inputs (per-gene association SUMMARY tables only -- no individual-level data), under RESULTS:
  <AD>/association_<tissue>_<trait>.tsv   per-gene TWAS z / p / predictor class, per tissue

Output: fig_assoc_combined.{pdf,png} in the shared figures directory.

Run:  GENE_ANNOT=<gencode.v19.genes.tsv> python analysis/figures/plot_assoc_combined.py

Input locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS); figures
go to results/figures/ by default (override with FLUX_FIGURES); the panel (c) Manhattan uses R
(ggplot2/ggrepel) when available and gene symbols from $GENE_ANNOT. Style: Arial, 300 DPI.
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.stats import norm
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils"))
import gene_labels as gl
import plot_manhattan as M

TIS = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
TRAIT = "CAD_aragam2022"
# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS = Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_cv"))
AD = f"{RESULTS}/horseshoe_alltargets/association"
# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
OUT = str(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))
os.makedirs(OUT, exist_ok=True)
ANNOT = gl.default_annot(os.environ.get("GENE_ANNOT"))
C_CIS, C_BOTH, C_TRANS = "#4C72B0", "#55A868", "#DD8452"
LAB = {"cis_only": "cis only", "cis_trans": "cis+trans", "trans_only": "trans only"}
COL = {"cis_only": C_CIS, "cis_trans": C_BOTH, "trans_only": C_TRANS}
NULL = (2 / np.pi) ** 0.5
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 9, "axes.linewidth": 0.8, "savefig.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False})

A = pd.concat([pd.read_csv(f"{AD}/association_{t}_{TRAIT}.tsv", sep="\t") for t in TIS], ignore_index=True)
sig = A[A["p_adj"] < 0.05].copy()
zcut = abs(norm.isf(sig["p"].max() / 2))

# ---- (c) Manhattan via the EXISTING R ggrepel renderer; title suppressed ----
M.ASSOC_DIR = Path(AD)
M.ASSOC_FMT = "association_{tissue}_" + TRAIT + ".tsv"
M.GENCODE = Path(ANNOT)
M.TITLE = "Transcriptome-wide association with coronary artery disease"
best = M.attach_positions(M.load_best_per_gene(), M.load_gene_positions())
n_sig_mapped = int((best["p_adj"] < M.FDR).sum())
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
m1 = (A["p_adj"] < 0.05) & (A["n_gwas"] == 1)                  # single-SNP (sign-flip-invariant)
axC.scatter(A.loc[m1, "sigma_g"], A.loc[m1, "z_twas"].abs(), s=30, facecolors="none",
            edgecolors="k", linewidths=0.8, label="single-SNP predictor")
axC.axhline(zcut, ls=":", lw=0.9, color="#c0392b")
r = np.corrcoef(az, 1.0 / A["sigma_g"].values)[0, 1]
axC.text(0.96, 0.96, f"$r(|Z_{{\\mathrm{{twas}}}}|,\\,1/\\sigma_g)={r:.2f}$", transform=axC.transAxes,
         ha="right", va="top", fontsize=7.5, color="#333")
axC.set_xlabel("$\\sigma_g$  (predicted-expression SD)")
axC.set_ylabel("$|Z_{\\mathrm{twas}}|$")
axC.set_title("(c) No $\\sigma_g$-driven inflation", loc="left", fontweight="bold")
axC.legend(frameon=False, fontsize=7, loc="lower right", handletextpad=0.2)

# ---- (d) embed the R Manhattan full-width, single panel title ----
axc = sf_bot.subplots(1, 1)
sf_bot.subplots_adjust(left=0.05, right=0.99, top=0.99, bottom=0.01)
axc.imshow(mpimg.imread(mh_png)); axc.axis("off")
axc.text(0.0, 1.0, "(d)", transform=axc.transAxes, fontsize=13, fontweight="bold",
         va="bottom", ha="left")

for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}/fig_assoc_combined.{ext}", bbox_inches="tight")
print("wrote", f"{OUT}/fig_assoc_combined.pdf  (manhattan R render ok={ok})")
print(f"sig gene-tissue={len(sig)} unique={sig.gene_id.nunique()} "
      f"trans-only(unique)={sig[sig.win_config=='trans_only'].gene_id.nunique()}")
