#!/usr/bin/env python3
"""Summary figure for the expression-prediction step: how well each gene's expression is
predicted from its cis SNPs and its trans (GRN-regulator) SNPs, and how much the trans channel
adds. One figure across all tissues, plus a per-tissue + pooled supplementary table.

Every expressed gene is a prediction target. Per gene, three horseshoe models are fit --
cis-only, trans-only and cis+trans -- and the best is selected by PSIS-LOO ELPD; a gene is
"predictable" when its selected model reaches cross-validated LOO-R^2 >= 0.01.

Inputs (per tissue; tissues are discovered from --models-dir unless --tissues is given):
  --models-dir  <tissue>.expr_model_stats.tsv      key->value summary (counts, per-channel
                                                   LOO-R^2 means) for the table and count panels
                <tissue>.expr_model_metrics.tsv.gz one row per gene (per-channel LOO-R^2, ELPD)
  --feat-dir    <tissue>.trans_features.tsv.gz     trans regulators per gene, for the
                                                   "number of regulator genes" axis

Output (--outdir):
  fig_alltargets_summary.{pdf,png}   five panels:
    (a) predictability funnel: genes modelled -> predictable -> trans-only -> predictable via trans
    (b) selected predictor per tissue (cis-only / trans-only / cis+trans), trans-using % on top
    (c) model evidence for trans: Delta-ELPD (cis+trans minus cis) histogram
    (d) trans-only predictable genes: LOO-R^2 vs number of regulator genes
    (e) who gains from trans: cis LOO-R^2 vs the LOO-R^2 gained by adding the trans channel
  table_alltargets_summary.tex       per-tissue + pooled supplementary table

Usage:
  plot_prediction_summary.py --models-dir DIR --feat-dir DIR --outdir DIR [--tissues T1 T2 ...]

Conventions: per-tissue breakdown only in figures (pooled totals in the funnel panel / table /
text), means not medians, Arial / 300 DPI / (a)(b) panel labels.
"""
import argparse
import glob
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8, "savefig.dpi": 300,
    "axes.spines.top": False, "axes.spines.right": False,
})

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--models-dir", required=True,
                help="dir with <tissue>.expr_model_stats.tsv and .expr_model_metrics.tsv.gz")
ap.add_argument("--feat-dir", required=True,
                help="dir with <tissue>.trans_features.tsv.gz (regulator-gene counts)")
ap.add_argument("--outdir", required=True, help="output dir for the figure and table")
ap.add_argument("--tissues", nargs="+", default=None,
                help="tissue labels in plotting order; default: discovered from --models-dir")
args = ap.parse_args()

MODELS = args.models_dir
FEAT = args.feat_dir
OUTDIR = args.outdir
os.makedirs(OUTDIR, exist_ok=True)
TISSUES = args.tissues or sorted(
    os.path.basename(f).split(".expr_model_stats.tsv")[0]
    for f in glob.glob(f"{MODELS}/*.expr_model_stats.tsv"))
if not TISSUES:
    raise SystemExit(f"no <tissue>.expr_model_stats.tsv found in {MODELS}")
TCOL = dict(zip(TISSUES, plt.cm.tab10(np.arange(len(TISSUES)))))
C_CIS, C_TRANS, C_BOTH = "#4C72B0", "#DD8452", "#55A868"


def load_stats(t):
    d = {}
    with open(f"{MODELS}/{t}.expr_model_stats.tsv") as fh:
        next(fh)
        for line in fh:
            k, v = line.rstrip("\n").split("\t")
            try:
                d[k] = float(v)
            except ValueError:
                d[k] = v
    return d

S = {t: load_stats(t) for t in TISSUES}

frames = []
for t in TISSUES:
    df = pd.read_csv(f"{MODELS}/{t}.expr_model_metrics.tsv.gz", sep="\t")
    df["tissue"] = t
    # --- schema bridge: new fit_expression_models output lacks delta_elpd ---
    df["delta_elpd"] = df["elpd_cistrans"] - df["elpd_cis"]
    # --- schema bridge: inject stats keys the new stats file does not carry ---
    # channel fit-counts (genes with a model fit in that channel) used by the
    # pooled-table weighted means; equal to the class counts in the new stats.
    S[t]["n_fit_cis"] = S[t]["n_cis_only"] + S[t]["n_both"]
    S[t]["n_fit_trans"] = S[t]["n_trans_only"] + S[t]["n_both"]
    S[t]["n_fit_cistrans"] = S[t]["n_both"]
    # mean best_loo_r2 over predictable genes (selected-model performance)
    S[t]["loo_r2_selected_predictable_mean"] = (
        df.loc[df["predictable"] == True, "best_loo_r2"].mean())
    frames.append(df)
M = pd.concat(frames, ignore_index=True)

# n_trans in the metrics counts trans SNPs (cis-SNPs of GRN ancestors). For the
# "# regulators" axis we want distinct ancestor GENES per target -> from trans_features.
reg = []
for t in TISSUES:
    tf = pd.read_csv(f"{FEAT}/{t}.trans_features.tsv.gz", sep="\t",
                     usecols=["gene_id", "source_gene_id"])
    g = tf.groupby("gene_id")["source_gene_id"].nunique().rename("n_reg_genes")
    g = g.reset_index(); g["tissue"] = t
    reg.append(g)
REG = pd.concat(reg, ignore_index=True)
M = M.merge(REG, on=["gene_id", "tissue"], how="left")

# fixed-seed shuffle so tissue colours overlay fairly in the scatter panels
M = M.sample(frac=1.0, random_state=0).reset_index(drop=True)
P = M[M["predictable"] == True].copy()            # predictable genes only

# ============================================================ figure
fig = plt.figure(figsize=(13.5, 8.4))
gs = fig.add_gridspec(2, 6, hspace=0.42, wspace=1.1)
axA = fig.add_subplot(gs[0, 0:2])   # (a) pooled funnel
axB = fig.add_subplot(gs[0, 2:4])   # (b) per-tissue composition
axC = fig.add_subplot(gs[0, 4:6])   # (c) trans-only predictable: LOO-R2 vs # regulator genes
axD = fig.add_subplot(gs[1, 0:3])   # (d) cis vs cis+trans LOO-R2 (does trans add)
axE = fig.add_subplot(gs[1, 3:6])   # (e) delta-ELPD model-evidence histogram


def pooled(metric):
    return int(sum(S[t][metric] for t in TISSUES))

# ---- (a) pooled funnel
g_all = pooled("n_genes")
g_pred = pooled("n_predictable")
g_to = pooled("n_trans_only")
g_to_pred = pooled("n_predictable_trans_only")
labels = ["Genes modeled",
          f"Predictable\n({100*g_pred/g_all:.0f}% of modeled)",
          f"Trans-only genes\n({100*g_to/g_all:.0f}%, no cis)",
          f"predictable via trans\n({100*g_to_pred/g_to:.0f}% of trans-only)"]
vals = [g_all, g_pred, g_to, g_to_pred]
y = np.arange(len(vals))[::-1]
axA.barh(y, vals, color="#9E9E9E", height=0.60)
for yi, v in zip(y, vals):
    axA.text(v + g_all * 0.015, yi, f"{v:,}", va="center", ha="left",
             fontsize=8.5, fontweight="bold")
axA.set_yticks(y); axA.set_yticklabels(labels, fontsize=7.8)
axA.set_xlim(0, g_all * 1.32)
axA.set_xlabel("Gene count (7 tissues pooled)")
axA.set_title("(a) Genome-wide gene predictability", loc="left", fontweight="bold")
axA.spines["left"].set_visible(False); axA.tick_params(left=False)

# ---- (b) per-tissue selected-model composition (predictable genes)
comp = (P.groupby(["tissue", "best_config"]).size()
          .unstack(fill_value=0).reindex(TISSUES))
x = np.arange(len(TISSUES))
b1 = comp["cis_only"].values
b2 = comp["trans_only"].values
b3 = comp["cis_trans"].values
axB.bar(x, b1, color=C_CIS, label="cis-only", width=0.72)
axB.bar(x, b2, bottom=b1, color=C_TRANS, label="trans-only", width=0.72)
axB.bar(x, b3, bottom=b1 + b2, color=C_BOTH, label="cis+trans", width=0.72)
tot = b1 + b2 + b3
for i in range(len(TISSUES)):
    axB.text(i, tot[i] + g_pred * 0.0016, f"{100*(b2[i]+b3[i])/tot[i]:.0f}%",
             ha="center", va="bottom", fontsize=7, color="#333333")
axB.set_xticks(x); axB.set_xticklabels(TISSUES)
axB.set_ylabel("Predictable genes")
axB.set_title("(b) Selected predictor per tissue", loc="left", fontweight="bold")
axB.legend(frameon=False, ncol=3, fontsize=7, loc="lower center",
           bbox_to_anchor=(0.5, -0.27), handlelength=1.0, columnspacing=1.0)


PT_COL = "#4C72B0"   # single colour for the per-gene scatter clouds

both = P[P["gene_class"] == "both"]
MS = 11   # scatter marker size (bigger = visible at panel scale)

# ---- (c) difference in ELPD (cis+trans - cis): direct model-evidence comparison.
# Selection is argmax ELPD, i.e. cis+trans wins iff Delta>0; Delta>4 = "clearly helps".
de = both["delta_elpd"]
XLO, XHI = -6, 30
bins = np.linspace(XLO, XHI, 49)
axC.hist(de[de <= 0].clip(XLO, XHI), bins=bins, color="#B0B0B0",
         edgecolor="white", linewidth=0.2, label="cis-only selected")
axC.hist(de[de > 0].clip(XLO, XHI), bins=bins, color=C_BOTH,
         edgecolor="white", linewidth=0.2, label="cis+trans selected")
axC.axvline(0, ls="--", lw=0.9, color="#333333")
axC.axvline(4, ls=":", lw=1.0, color="#A24C3D")
axC.set_xlim(XLO, XHI)
axC.set_xlabel(r"$\Delta$ELPD  (cis+trans $-$ cis)")
axC.set_ylabel("Both-channel genes")
axC.set_title("(c) Model evidence for trans", loc="left", fontweight="bold")
axC.legend(frameon=False, fontsize=7, loc="upper right", handlelength=1.0)
axC.text(0.96, 0.68, f"{100*(de>0).mean():.0f}% $\\Delta>0$\n(trans selected)",
         transform=axC.transAxes, va="top", ha="right", fontsize=7.5)
axC.text(4, axC.get_ylim()[1] * 0.55, " $\\Delta$=4", fontsize=6.8,
         color="#A24C3D", va="center", ha="left")

# ---- (d) trans-only genes: LOO-R2 vs number of regulator GENES (GRN ancestors)
to = P[P["gene_class"] == "trans_only"]
axD.scatter(to["n_reg_genes"], to["loo_r2_trans"], s=MS, alpha=0.30,
            c=PT_COL, rasterized=True, linewidths=0)
axD.axhline(0.01, ls=":", lw=0.8, color="#999999")
axD.set_xlim(1, to["n_reg_genes"].max() * 1.02)   # >=1 regulator (0 parents is meaningless)
axD.set_xlabel("# trans regulator genes (GRN ancestors)")
axD.set_ylabel("trans-only LOO-$R^2$")
axD.set_title("(d) Trans-only predictable genes", loc="left", fontweight="bold")
axD.text(0.96, 0.96, f"n = {len(to):,} genes (no cis signal)",
         transform=axD.transAxes, va="top", ha="right", fontsize=7.5)

# ---- (e) who gains from trans: cis LOO-R2 vs delta LOO-R2 (both-channel genes)
delta = both["loo_r2_cistrans"] - both["loo_r2_cis"]
axE.scatter(both["loo_r2_cis"], delta, s=MS, alpha=0.30,
            c=PT_COL, rasterized=True, linewidths=0)
axE.axhline(0, ls="--", lw=0.8, color="#444444")
axE.set_ylim(0, 0.2)                          # zoom to the informative gain range
axE.set_xlim(0, both["loo_r2_cis"].max() * 1.02)
axE.set_xlabel("cis-only LOO-$R^2$")
axE.set_ylabel(r"$\Delta$ LOO-$R^2$ (cis+trans $-$ cis)")
axE.set_title("(e) Who gains from trans?", loc="left", fontweight="bold")

for ext in ("pdf", "png"):
    fig.savefig(f"{OUTDIR}/fig_alltargets_summary.{ext}", bbox_inches="tight")
print(f"wrote {OUTDIR}/fig_alltargets_summary.pdf / .png")

# ============================================================ LaTeX table
def wmean(mkey, nkey):
    num = sum(S[t][mkey] * S[t][nkey] for t in TISSUES)
    den = sum(S[t][nkey] for t in TISSUES)
    return num / den

rows = []
for t in TISSUES:
    s = S[t]
    su = 100 * (s["n_selected_trans_only"] + s["n_selected_cis_trans"]) / s["n_predictable"]
    rows.append((t, int(s["n_genes"]), int(s["n_predictable"]),
                 100 * s["n_predictable"] / s["n_genes"],
                 int(s["n_predictable_trans_only"]), su,
                 s["loo_r2_cis_mean"], s["loo_r2_trans_mean"],
                 s["loo_r2_cistrans_mean"], s["loo_r2_selected_predictable_mean"]))
pg, pp, ppto = pooled("n_genes"), pooled("n_predictable"), pooled("n_predictable_trans_only")
psu = 100 * (pooled("n_selected_trans_only") + pooled("n_selected_cis_trans")) / pp

tex = [r"\begin{table}[ht]", r"\centering",
       r"\caption{Genome-wide (all-targets) horseshoe expression models across seven "
       r"STARNET tissues. Every expressed gene is modelled as a prediction target from its "
       r"cis SNPs and its trans (GRN-parent) regulators; the per-gene predictor "
       r"(cis-only, trans-only, or cis+trans) is selected by PSIS-LOO ELPD. A gene is "
       r"predictable when its selected model attains cross-validated LOO-$R^2\geq0.01$. "
       r"Per-channel LOO-$R^2$ means are over genes fit in that channel. All summaries are means.}",
       r"\label{tab:alltargets_summary}", r"\small",
       r"\begin{tabular}{lrrrrrcccc}", r"\toprule",
       r" & & & & Trans-only & Trans- & \multicolumn{4}{c}{Mean LOO-$R^2$} \\",
       r"\cmidrule(lr){7-10}",
       r"Tissue & Genes & Predict. & \% pred. & predict. & using \% & cis & trans & cis+trans & selected \\",
       r"\midrule"]
for (t, g, p, pc, pto, su, rc, rt, rct, rs) in rows:
    tex.append(f"{t} & {g:,} & {p:,} & {pc:.1f} & {pto:,} & {su:.0f} & "
               f"{rc:.3f} & {rt:.3f} & {rct:.3f} & {rs:.3f} \\\\")
tex += [r"\midrule",
        f"\\textbf{{Pooled}} & {pg:,} & {pp:,} & {100*pp/pg:.1f} & {ppto:,} & {psu:.0f} & "
        f"{wmean('loo_r2_cis_mean','n_fit_cis'):.3f} & "
        f"{wmean('loo_r2_trans_mean','n_fit_trans'):.3f} & "
        f"{wmean('loo_r2_cistrans_mean','n_fit_cistrans'):.3f} & "
        f"{wmean('loo_r2_selected_predictable_mean','n_predictable'):.3f} \\\\",
        r"\bottomrule", r"\end{tabular}", r"\end{table}"]
tabpath = f"{OUTDIR}/table_alltargets_summary.tex"
with open(tabpath, "w") as fh:
    fh.write("\n".join(tex) + "\n")
print(f"wrote {tabpath}")
print("\nPooled: {:,} genes | {:,} predictable ({:.1f}%) | {:,} trans-only-predictable "
      "| {:.0f}% trans-using".format(pg, pp, 100 * pp / pg, ppto, psu))
