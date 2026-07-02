#!/usr/bin/env python3
"""Functional-relevance axis figure: are flux-network gene classes functionally
cardiovascular? Forest plot of enrichment odds ratios (95% CI) vs the tested-gene
background, for (a) a curated cardiovascular gene universe (MGI mouse-KO phenotype +
GO-BP + Reactome + KEGG + Jensen DISEASES) and (b) mouse-knockout cardiovascular
phenotypes alone (functionally orthogonal to genetics). Headline: trans_only alone
is weak, but CONVERGENCE upgrades it -- convergent trans_only genes are CV-enriched.

Inputs are per-gene SUMMARY tables under RESULTS (associations and flux node roles) plus the
curated gene-set libraries shipped in data/genesets/. Output: fig_functional_axis.{pdf,png} in the
shared figures directory.

Run:  GENE_ANNOT=<gencode.v19.genes.tsv> python analysis/figures/plot_functional_axis.py

Input locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS); figures
go to results/figures/ by default (override with FLUX_FIGURES); gene symbols come from the
annotation in $GENE_ANNOT. Style: Arial, 300 DPI.
"""
import sys, re, os, numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "utils")); import gene_labels as gl
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
                     "font.size":9,"axes.linewidth":0.8,"savefig.dpi":300,"pdf.fonttype":42})
# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS=Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_cv"))
AT=f"{RESULTS}/horseshoe_alltargets"; TIS=["AOR","Blood","LIV","MAM","SF","SKLM","VAF"]
# Curated gene-set libraries (.gmt) ship in the repo under data/genesets/.
TRAIT="CAD_aragam2022"; GS=str(Path(__file__).resolve().parents[2] / "data" / "genesets")
# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
OUT=str(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))
os.makedirs(OUT, exist_ok=True)
base=lambda g:str(g).split(".")[0]
MAP=gl.load_label_map(gl.default_annot(os.environ.get("GENE_ANNOT"))); sym=lambda g:gl.resolve(g,MAP)
CV=re.compile(r"cardiovascular|cardiac|heart|myocard|cardiomyo|vascul|vasocon|vasodil|"
  r"angiogen|endotheli|blood vessel|artery|arterial|aorta|aortic|atheroscler|"
  r"smooth muscle|blood pressure|hypertens|coagul|platelet|thrombo|hemostasis|fibrin|"
  r"lipid|lipoprotein|cholesterol|choleste|triglycer|sterol|foam cell|"
  r"coronary|ischemi|infarct|elastic fiber|extracellular matrix",re.I)
def load_gmt(p,pf):
    d={}
    for ln in open(p):
        f=ln.rstrip("\n").split("\t")
        if len(f)>=3: d[pf+":"+f[0]]=set(g.upper() for g in f[2:] if g)
    return d
LIBS={"MGI":load_gmt(f"{GS}/MGI_Mammalian_Phenotype_Level_4_2021.gmt","MGI"),
      "GO-BP":load_gmt(f"{GS}/GO_Biological_Process_2023.gmt","GOBP"),
      "REAC":load_gmt(f"{GS}/Reactome_2022.gmt","REAC"),"KEGG":load_gmt(f"{GS}/KEGG_2021_Human.gmt","KEGG"),
      "Jensen":load_gmt(f"{GS}/Jensen_DISEASES.gmt","DIS")}

A=pd.concat([pd.read_csv(f"{AT}/association/association_{t}_{TRAIT}.tsv",sep="\t") for t in TIS],ignore_index=True)
A["b"]=A.gene_id.map(base); A=A.sort_values("p_adj"); best=A.drop_duplicates("b",keep="first")
tested=set(A.b.unique()); sig=best[best.p_adj<0.05]
S=pd.read_csv(f"{AT}/functional/functional_relevance_summary.tsv",sep="\t")
conv=set(S.loc[S.convergent,"tbase"])
cisonly=set(sig.loc[sig.win_config=="cis_only","b"]); cistrans=set(sig.loc[sig.win_config=="cis_trans","b"])
transonly=set(sig.loc[sig.win_config=="trans_only","b"])
tested_sym=set(sym(g) for g in tested)
CVgenes=set().union(*[gs for d in LIBS.values() for t,gs in d.items() if CV.search(t)])&tested_sym
MGIcv=set().union(*[gs for t,gs in LIBS["MGI"].items() if CV.search(t)])&tested_sym

CLASSES=[("cis_only",cisonly,"#7f7f7f"),("cis+trans",cistrans,"#4C72B0"),
         ("trans_only",transonly,"#dd8452"),("convergent\n(all cores)",conv,"#2a8f6b"),
         ("convergent\n& trans_only",conv&transonly,"#2a8f6b"),
         ("convergent\n& cis+trans",conv&cistrans,"#2a8f6b")]
def OR_ci(geneset,universe):
    qs=set(sym(g) for g in geneset)&tested_sym; n=len(qs)
    a=len(qs&universe); b=n-a; rest=tested_sym-qs; c=len(rest&universe); d=len(rest)-c
    _,p=stats.fisher_exact([[a,b],[c,d]],alternative="greater")
    aa,bb,cc,dd=[x+0.5 for x in (a,b,c,d)]; orr=(aa*dd)/(bb*cc)
    se=np.sqrt(1/aa+1/bb+1/cc+1/dd); lo,hi=np.exp(np.log(orr)-1.96*se),np.exp(np.log(orr)+1.96*se)
    return orr,lo,hi,p,n,100*a/n
pfmt=lambda p:(f"p={p:.0e}" if p<1e-3 else f"p={p:.2g}") if p<.05 else f"p={p:.2g} (n.s.)"

fig,(axA,axB)=plt.subplots(1,2,figsize=(11,4.6),sharey=True)
y=np.arange(len(CLASSES))[::-1]
for ax,uni,title,bgset in [(axA,CVgenes,"(a) Cardiovascular gene universe","CVgenes"),
                           (axB,MGIcv,"(b) Mouse-knockout cardiovascular phenotype","MGIcv")]:
    bgfrac=100*len(uni)/len(tested_sym)
    for yi,(lab,gs,col) in zip(y,CLASSES):
        orr,lo,hi,p,n,pct=OR_ci(gs,uni)
        ax.plot([lo,hi],[yi,yi],color=col,lw=1.6,zorder=2)
        ax.scatter([orr],[yi],s=58,color=col,edgecolors="white",linewidths=0.6,zorder=3)
        ax.text(hi+0.04,yi,pfmt(p),va="center",ha="left",fontsize=7.2,color=col)
    ax.axvline(1,ls="--",lw=1,color="#888",zorder=1)
    ax.text(0.99,np.mean(y),"OR = 1  (no enrichment)",rotation=90,ha="right",va="center",fontsize=7,color="#888")
    ax.set_xlabel("enrichment odds ratio  (95% CI)")
    ax.set_title(title,loc="left",fontweight="bold",fontsize=9.6)
    ax.set_xlim(0.7,None); ax.spines[["top","right"]].set_visible(False); ax.grid(axis="x",color="0.93",lw=.5)
axA.set_yticks(y); axA.set_yticklabels([c[0] for c in CLASSES],fontsize=8)
fig.tight_layout()
for ext in ("pdf","png"): fig.savefig(f"{OUT}/fig_functional_axis.{ext}",bbox_inches="tight")
print("classes (CV-universe): "+" | ".join(f"{c[0].replace(chr(10),' ')} OR={OR_ci(c[1],CVgenes)[0]:.2f} {pfmt(OR_ci(c[1],CVgenes)[3])}" for c in CLASSES))
print(f"wrote {OUT}/fig_functional_axis.pdf")
