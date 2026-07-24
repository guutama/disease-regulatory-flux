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

Run:  GENE_ANNOT=<gencode.v19.genes.tsv> python analysis/figure_scripts/plot_functional_axis.py

Input locations resolve under RESULTS (default results_findr_py/, override with FLUX_RESULTS); figures
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
# Input result locations resolve under RESULTS (default results_findr_py/, override with FLUX_RESULTS).
# association and functional tables live directly under results_findr_py/.
RESULTS=Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_findr_py"))
AT=f"{RESULTS}"; TIS=["AOR","Blood","LIV","MAM","SF","SKLM","VAF"]
# Curated gene-set libraries (.gmt) ship in the repo under data/genesets/.
TRAIT="CAD_aragam2022"; GS=str(Path(__file__).resolve().parents[2] / "data" / "genesets")
OUT=str(os.environ.get("FLUX_FIGURES", RESULTS / "figures"))
ANNOT=os.environ.get("GENE_ANNOT",
    "/cluster/projects/nn1015k/datasets/references/gencode/GRCh37/gencode.v19.genes.tsv")
os.makedirs(OUT, exist_ok=True)
base=lambda g:str(g).split(".")[0]
MAP=gl.load_label_map(ANNOT); sym=lambda g:gl.resolve(g,MAP)
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
# rename the selected-config column config -> win_config for the code below
if "win_config" not in A.columns and "config" in A.columns:
    A["win_config"]=A["config"]
A["b"]=A.gene_id.map(base); A=A.sort_values("p_adj"); best=A.drop_duplicates("b",keep="first")
tested=set(A.b.unique()); sig=best[best.p_adj<0.05]
S=pd.read_csv(f"{AT}/functional/functional_relevance_summary.tsv",sep="\t")
conv=set(S.loc[S.convergent,"tbase"])
cisonly=set(sig.loc[sig.win_config=="cis_only","b"])
transusing=set(sig.loc[sig.win_config.isin(["cis_trans","trans_only"]),"b"])
# non-significant contributors: regulators that feed trans signal but are not themselves disease genes
_tru=lambda x:str(x).strip().lower() in ("true","1"); noncontrib=set()
for t in TIS:
    nd=pd.read_csv(f"{AT}/flux/nodes/flux_nodes_{t}_{TRAIT}.tsv",sep="\t")
    nd["b"]=nd.gene_id.map(base); noncontrib|=set(nd.loc[~nd.disease_sig.map(_tru),"b"])
tested_sym=set(sym(g) for g in tested)
CVgenes=set().union(*[gs for d in LIBS.values() for t,gs in d.items() if CV.search(t)])&tested_sym
MGIcv=set().union(*[gs for t,gs in LIBS["MGI"].items() if CV.search(t)])&tested_sym

CLASSES=[("cis_only",cisonly,"#7f7f7f"),("trans_using",transusing,"#4C72B0"),
         ("convergent",conv,"#2a8f6b"),("non-sig.\ncontributor",noncontrib,"#9467bd")]
def OR_ci(geneset,universe):
    qs=set(sym(g) for g in geneset)&tested_sym; n=len(qs)
    a=len(qs&universe); b=n-a; rest=tested_sym-qs; c=len(rest&universe); d=len(rest)-c
    _,p=stats.fisher_exact([[a,b],[c,d]],alternative="greater")
    aa,bb,cc,dd=[x+0.5 for x in (a,b,c,d)]; orr=(aa*dd)/(bb*cc)
    se=np.sqrt(1/aa+1/bb+1/cc+1/dd); lo,hi=np.exp(np.log(orr)-1.96*se),np.exp(np.log(orr)+1.96*se)
    return orr,lo,hi,p,n,100*a/n
pfmt=lambda p:(f"p={p:.0e}" if p<1e-3 else f"p={p:.2g}") if p<.05 else f"p={p:.2g} (n.s.)"

# cardiovascular terms (for panel c): term -> gene set, and a cleaned display name
CVTERMS={t:gs for d in LIBS.values() for t,gs in d.items() if CV.search(t)}
def clean_term(t):
    t=re.sub(r'^[A-Za-z\-]+:','',t)                                   # drop library prefix
    t=re.sub(r'\s*\(?(GO|MP|R-HSA|WP|HP|DOID):?[0-9A-Za-z\-]+\)?\s*$','',t)  # drop trailing ontology id
    return t.strip()
def top_terms(geneset,topn=6,kmin=3,Kmin=5,Kmax=1500):
    """Cardiovascular terms most over-represented in a class, by fold enrichment (observed/expected)."""
    qs=set(sym(g) for g in geneset)&tested_sym; n=len(qs); N=len(tested_sym)
    rows=[]
    for t,gs in CVTERMS.items():
        gsT=gs&tested_sym; K=len(gsT)
        if K<Kmin or K>Kmax: continue
        k=len(gsT&qs)
        if k<kmin: continue
        rows.append((clean_term(t),k,K,k/(n*K/N)))
    return pd.DataFrame(rows,columns=["term","k","K","fold"]).sort_values("fold",ascending=False).head(topn)

# random size-matched null: a random set of tested genes should give OR ~ 1 (no enrichment).
TESTED_LIST=sorted(tested_sym); RNG=np.random.default_rng(0)
def null_band(n,universe,B=1000):
    """95% envelope of the enrichment OR for random n-gene sets drawn from the tested-gene background."""
    cvmask=np.fromiter((g in universe for g in TESTED_LIST),bool,len(TESTED_LIST))
    N=len(TESTED_LIST); Ncv=int(cvmask.sum()); ors=np.empty(B)
    for b in range(B):
        a=int(cvmask[RNG.choice(N,size=n,replace=False)].sum())
        aa,bb,cc,dd=a+.5,n-a+.5,Ncv-a+.5,(N-n)-(Ncv-a)+.5; ors[b]=(aa*dd)/(bb*cc)
    return np.percentile(ors,[2.5,97.5])

import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
fig=plt.figure(figsize=(11,8.6))
gsp=gridspec.GridSpec(2,2,height_ratios=[1,2.5],hspace=0.45,wspace=0.22)
axA=fig.add_subplot(gsp[0,0]); axB=fig.add_subplot(gsp[0,1],sharey=axA); axC=fig.add_subplot(gsp[1,:])
y=np.arange(len(CLASSES))[::-1]
for ax,uni,title in [(axA,CVgenes,"(a) Cardiovascular gene universe"),
                     (axB,MGIcv,"(b) Mouse-knockout cardiovascular phenotype")]:
    for yi,(lab,gs,col) in zip(y,CLASSES):
        orr,lo,hi,p,n,pct=OR_ci(gs,uni)
        nlo,nhi=null_band(len(set(sym(g) for g in gs)&tested_sym),uni)
        ax.barh(yi,nhi-nlo,left=nlo,height=0.52,color="#e6e6e6",zorder=0)   # random size-matched null (95%)
        ax.plot([lo,hi],[yi,yi],color=col,lw=1.6,zorder=2)
        ax.scatter([orr],[yi],s=58,color=col,edgecolors="white",linewidths=0.6,zorder=3)
        ax.text(hi+0.04,yi,pfmt(p),va="center",ha="left",fontsize=7.2,color=col)
    ax.axvline(1,ls="--",lw=1,color="#888",zorder=1)
    ax.set_xlabel("enrichment odds ratio  (95% CI)")
    ax.set_title(title,loc="left",fontweight="bold",fontsize=9.6)
    ax.set_xlim(0.7,None); ax.spines[["top","right"]].set_visible(False); ax.grid(axis="x",color="0.93",lw=.5)
axA.set_yticks(y); axA.set_yticklabels([c[0] for c in CLASSES],fontsize=8)
axB.tick_params(labelleft=False)
axA.legend(handles=[Patch(color="#e6e6e6",label="random size-matched null (95%)")],
           frameon=False,fontsize=6.8,loc="upper center",bbox_to_anchor=(0.5,-0.28),ncol=1)

# (c) most-enriched cardiovascular terms per class, by fold enrichment
rows=[]; yp=0.0
for lab,gs,col in CLASSES:
    for _,r in top_terms(gs).iterrows():
        rows.append((yp,r.fold,col,r.term,int(r.k),int(r.K))); yp-=1.0
    yp-=0.9                                                            # gap between class blocks
for ypi,fold,col,term,k,K in rows:
    axC.barh(ypi,fold,color=col,edgecolor="white",height=0.74,zorder=3)
    axC.text(fold+0.06,ypi,f"{k}/{K}",va="center",ha="left",fontsize=6.2,color="#888")
axC.axvline(1,ls="--",lw=1,color="#888",zorder=1)
axC.set_yticks([r[0] for r in rows])
axC.set_yticklabels([(r[3][:46]+"…") if len(r[3])>46 else r[3] for r in rows],fontsize=6.9)
axC.set_xlabel("fold enrichment  (observed / expected)"); axC.set_xlim(0,None)
axC.set_title("(c) Most-enriched cardiovascular terms per gene class",
              loc="left",fontweight="bold",fontsize=9.6)
axC.spines[["top","right","left"]].set_visible(False); axC.tick_params(axis="y",length=0)
axC.grid(axis="x",color="0.93",lw=.5)
axC.legend(handles=[Patch(color=c[2],label=c[0]) for c in CLASSES],frameon=False,fontsize=8,loc="lower right")
for ext in ("pdf","png"): fig.savefig(f"{OUT}/fig_functional_axis.{ext}",bbox_inches="tight")
print("classes (CV-universe): "+" | ".join(f"{c[0]} OR={OR_ci(c[1],CVgenes)[0]:.2f} {pfmt(OR_ci(c[1],CVgenes)[3])}" for c in CLASSES))
print(f"wrote {OUT}/fig_functional_axis.pdf")
