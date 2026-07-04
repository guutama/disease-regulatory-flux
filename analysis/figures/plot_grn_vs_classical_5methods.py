#!/usr/bin/env python3
"""Five-method prediction comparison figure (``fig_grn_vs_classical_5methods``).

Two panels comparing three GRN-TWAS expression priors (Bayesian ridge, BSLMM,
regularised horseshoe) against two classical cis-only predictors (Elastic-net /
PrediXcan, GBLUP), per tissue:
  (a) number of predictable genes (held-out R2 >= 0.01) per method;
  (b) mean held-out R2 on the cis-eGenes shared by all methods (a like-for-like set).

Input: a precomputed per-tissue counts table (comparison5_counts.tsv) under RESULTS.
Output: fig_grn_vs_classical_5methods.{pdf,png} in the shared figures directory.

Input locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS); figures
go to results/figures/ by default (override with FLUX_FIGURES). Style: Arial, 300 DPI.
"""
import os
from pathlib import Path

import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
                     "font.size":9,"axes.spines.top":False,"axes.spines.right":False,"axes.axisbelow":True})
# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS = Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_cv"))
# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
OUT = str(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))
os.makedirs(OUT, exist_ok=True)
df=pd.read_csv(f"{RESULTS}/horseshoe_alltargets/figures/comparison5_counts.tsv",sep="\t")
T=df.tissue.tolist(); x=np.arange(len(T)); w=0.16
# order: enet, GBLUP, Bayesian ridge, BSLMM, horseshoe
C_EN="#3b6ea5"; C_GB="#9ec3e0"; C_BR="#80cdc1"; C_BS="#2e8b57"; C_HS="#00441b"
pos=[-2*w,-w,0,w,2*w]
SPEC=[("enet","pred","r2_enet",C_EN,"Elastic-net (PrediXcan)","Elastic-net"),
      ("gblup","pred","r2_gblup",C_GB,"GBLUP (ridge)","GBLUP"),
      ("br","pred","r2_br",C_BR,"GRN-TWAS: Bayesian ridge","Bayesian ridge"),
      ("bs","pred","r2_bs",C_BS,"GRN-TWAS: BSLMM","BSLMM"),
      ("hs","pred","r2_hs",C_HS,"GRN-TWAS: horseshoe","horseshoe")]
fig,(axA,axB)=plt.subplots(1,2,figsize=(12.5,4.3),gridspec_kw={"wspace":0.24})
for p,(pre,_,r2,c,_,_) in zip(pos,SPEC):
    axA.bar(x+p, df[f"{pre}_pred"], w, color=c, edgecolor="white", linewidth=.4)
    axB.bar(x+p, df[r2], w, color=c, edgecolor="white", linewidth=.4)
axA.set_xticks(x); axA.set_xticklabels(T); axB.set_xticks(x); axB.set_xticklabels(T)
axA.set_ylabel("predictable genes  ($R^2\\geq0.01$)"); axA.set_ylim(0, df.hs_pred.max()*1.5)
axB.set_ylabel("mean held-out $R^2$  (shared cis-eGenes)"); axB.set_ylim(0, df.r2_hs.max()*1.34)
axA.set_title("(a) Predictable-gene count", loc="left", fontsize=10.5, fontweight="bold")
axB.set_title("(b) Accuracy on shared cis-eGenes", loc="left", fontsize=10.5, fontweight="bold")
for ax in (axA,axB): ax.grid(axis="y", color="0.9", linewidth=.6)
axA.legend(handles=[Patch(fc=c,label=l) for *_,c,l,_ in SPEC], fontsize=7, frameon=False,
           loc="upper left", ncol=2, columnspacing=1.1, handlelength=1.4)
axB.legend(handles=[Patch(fc=c,label=s) for *_,c,_,s in SPEC], fontsize=7, frameon=False,
           loc="upper left", ncol=2, columnspacing=1.1, handlelength=1.4)
for ext in ("png","pdf"):
    fig.savefig(f"{OUT}/fig_grn_vs_classical_5methods.{ext}", dpi=300, bbox_inches="tight")
print(f"wrote {OUT}/fig_grn_vs_classical_5methods.pdf")
