#!/usr/bin/env python3
"""Main-readouts figure of the disease-regulatory-flux manuscript (``fig_main_readouts``).

This script reproduces the six-panel main-text figure that reads a gene's disease role off the
reconstructed regulatory-flux network along four axes:
  (a) Direct evidence   -- the gene's own genetic signal |z_TWAS| (is it itself a disease gene);
  (b) Involvement       -- how much it receives (in-degree) vs broadcasts (out-degree), by class;
  (c) Key drivers       -- non-significant broadcasters: own signal vs out-flux broadcast breadth;
  (d) Core master regs  -- significant broadcasters: own signal vs out-flux broadcast breadth;
  (e) Coordinated OUT   -- do broadcasters drive their targets coherently (out-coherence > 0.7)?
  (f) Convergence IN    -- are disease genes regulated coherently (in-coherence > 0.7)?
Panels (e) and (f) share the coherence axis to show the asymmetry (broadcasting is incoherent,
regulation of disease genes is convergent); each shows a sign-flip permutation null with the
observed value drawn on top.

Inputs (per-gene / per-edge SUMMARY tables only -- no individual-level data):
  <FLUX>/flux_nodes_<tissue>_<trait>.tsv   per-gene node roles, flux and coherence terms
  <FLUX>/flux_edges_<tissue>_<trait>.tsv   regulator -> target flux edges (for the coherence panels)

Outputs (written to <OUTDIR>):
  fig_main_readouts.{pdf,png}              the six-panel figure
  <FLUX>/perm_null_{in,out}_<trait>.npz    cached sign-flip null (reused across renders)

Run:  GENE_ANNOT=<gencode.v19.genes.tsv> python analysis/figure_scripts/plot_main_readouts.py

Input result locations resolve under RESULTS -- by default the repo's own ``results_findr_py/``
directory, overridable with the ``FLUX_RESULTS`` environment variable. It reads the flux node/edge
tables under ``<RESULTS>/flux/nodes`` and writes figures to ``<RESULTS>/figures`` (override with
``FLUX_FIGURES``). Gene symbols are resolved from ``$GENE_ANNOT``. Set ``RECOMPUTE_NULL=1`` (or
delete the .npz) to force a fresh permutation. Style: Arial, 300 DPI, (a)(b) panel labels.
"""
import pandas as pd, numpy as np, sys, os
from pathlib import Path
_REPO = Path(__file__).resolve().parents[2]                  # analysis/figure_scripts/ -> repo root
sys.path.insert(0, str(_REPO / "utils")); import gene_labels as gl
from adjustText import adjust_text
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
                     "font.size":9,"axes.linewidth":0.8,"savefig.dpi":300,"pdf.fonttype":42})
# Result files this figure reads/writes, resolved under RESULTS (default: the repo's
# results_findr_py/; override the root with the FLUX_RESULTS environment variable). Reads the
# flux node/edge tables under <RESULTS>/flux/nodes, writes figures to <RESULTS>/figures.
RESULTS = Path(os.environ.get("FLUX_RESULTS", _REPO / "results_findr_py"))
FL = RESULTS / "flux" / "nodes"; AT = FL
OUTDIR = Path(os.environ.get("FLUX_FIGURES", RESULTS / "figures"))
ANNOT = os.environ.get("GENE_ANNOT",
    "/cluster/projects/nn1015k/datasets/references/gencode/GRCh37/gencode.v19.genes.tsv")
os.makedirs(OUTDIR, exist_ok=True)
TIS=["AOR","Blood","LIV","MAM","SF","SKLM","VAF"]
base=lambda g:str(g).split(".")[0]; THR=0.7
N=pd.concat([pd.read_csv(f"{FL}/flux_nodes_{t}_CAD_aragam2022.tsv",sep="\t") for t in TIS],ignore_index=True)
# out_degree is the broadcast breadth (# distinct disease targets fed); alias it to
# n_disease_targets, the name the panels below use.
if "n_disease_targets" not in N.columns:
    N["n_disease_targets"]=N["out_degree"]
for c in ["z_twas","out_flux","net_out","in_flux","net_in","in_degree","out_degree","n_disease_targets"]:
    N[c]=pd.to_numeric(N[c],errors="coerce")
N["az"]=N.z_twas.abs().fillna(0); N["sig"]=N.disease_sig==True
N["coh_out"]=(N.net_out.abs()/N.out_flux.replace(0,np.nan))
N["coh_in"] =(N.net_in.abs()/N.in_flux.replace(0,np.nan))
MAP=gl.load_label_map(ANNOT)
N["sym"]=N.gene_id.map(lambda g:gl.resolve(g,MAP))
N["pmr"]=N["is_pmr"].astype(str).str.lower().isin(["true","1","1.0"])
thr=N.loc[N.sig,"az"].min()
# Shared own-signal |z_TWAS| ceiling for panels a, c, d, so the z-score axis is on one scale
# across all three. Set by the broadcasters shown in c/d; a handful of very
# high-|z| disease genes beyond it are clipped into panel a's top bin (noted on the axis).
_bmask=N.n_disease_targets.fillna(0)>10
ZMAX=float(np.ceil(N.loc[_bmask,"az"].max()))+1
BOX=dict(boxstyle="round,pad=0.3",fc="white",ec="#cccccc",alpha=0.92)

fig,ax=plt.subplots(2,3,figsize=(16.5,9)); (A,B,C),(D,E,F)=ax

# (a) DIRECT EVIDENCE -- own-signal distribution on the shared |z_TWAS| scale (0..ZMAX);
# the few genes with |z|>ZMAX fall outside the axis range and are not shown.
A.hist(N.loc[N.sig,"az"],bins=np.linspace(0,ZMAX,45),color="#c1121f",alpha=.85,label=f"disease gene  ({int(N.sig.sum()):,})")
A.hist(N.loc[~N.sig,"az"],bins=np.linspace(0,ZMAX,45),color="#cccccc",alpha=.85,label=f"not significant  ({int((~N.sig).sum()):,})")
A.axvline(thr,ls="--",lw=1,color="#333")
A.set_ylim(0,A.get_ylim()[1]*1.12)
A.annotate(f"FDR 0.05\n|z| = {thr:.1f}",xy=(thr,A.get_ylim()[1]*.5),xytext=(thr+1.6,A.get_ylim()[1]*.62),
           fontsize=7.5,color="#333",ha="center",bbox=BOX,arrowprops=dict(arrowstyle="->",color="#333",lw=.8))
A.set_xlim(0,ZMAX); A.set_xlabel("own genetic signal  |z$_{\\mathrm{TWAS}}$|  (z-score)"); A.set_ylabel("genes"); A.legend(frameon=False,fontsize=8,loc="upper right")
A.set_title("(a) Direct evidence — is the gene a disease gene?",loc="left",fontweight="bold",fontsize=9.7)

# (b) INVOLVEMENT = receiving (in-degree) vs broadcasting (out-degree), per gene class -- diverging bars
N["bcast"]=N.n_disease_targets>10
def _cls(r): return ("core master reg" if r.sig and r.bcast else "core (sink)" if r.sig else "key driver" if r.bcast else "minor regulator")
N["cls"]=N.apply(_cls,axis=1)
classes=["minor regulator","key driver","core (sink)","core master reg"]; y=np.arange(len(classes))
inm =[N[N.cls==q].in_degree.fillna(0).mean()  for q in classes]
outm=[N[N.cls==q].out_degree.fillna(0).mean() for q in classes]
insd =[N[N.cls==q].in_degree.fillna(0).std(ddof=0)  for q in classes]   # SD across gene-tissue instances
outsd=[N[N.cls==q].out_degree.fillna(0).std(ddof=0) for q in classes]
B.barh(y,[-v for v in inm],xerr=insd,color="#2166ac",edgecolor="white",height=.7,
       error_kw=dict(ecolor="#0b355e",elinewidth=.9,capsize=2.5),label="receives from  (# regulators)")
B.barh(y,outm,xerr=outsd,color="#e08214",edgecolor="white",height=.7,
       error_kw=dict(ecolor="#8a4b06",elinewidth=.9,capsize=2.5),label="broadcasts to  (# disease genes)")
for yi,(iv,ivs,ov,ovs) in enumerate(zip(inm,insd,outm,outsd)):
    if iv>0.3: B.text(-iv-ivs-0.8,yi,f"{iv:.1f}$\\pm${ivs:.1f}",ha="right",va="center",fontsize=6.8)
    if ov>0.3: B.text(ov+ovs+0.8,yi,f"{ov:.1f}$\\pm${ovs:.1f}",ha="left",va="center",fontsize=6.8)
B.axvline(0,color="#333",lw=.9)
B.set_yticks(y); B.set_yticklabels([f"{c}\n(n={int((N.cls==c).sum())})" for c in classes],fontsize=8)
B.set_xlabel("$\\leftarrow$ receives from        broadcasts to $\\rightarrow$   (mean $\\pm$ SD # links)")
mx=max(max(i+s for i,s in zip(inm,insd)),max(o+s for o,s in zip(outm,outsd)))*1.32; B.set_xlim(-mx,mx)
from matplotlib.ticker import FuncFormatter
B.xaxis.set_major_formatter(FuncFormatter(lambda v,_: f"{abs(v):g}"))   # show |value| (both sides positive)
B.legend(frameon=False,fontsize=7.3,loc="lower left")
B.set_title("(b) Involvement — receiving vs broadcasting",loc="left",fontweight="bold",fontsize=9.7)

# (c) COORDINATED OUT  vs  (d) CONVERGENCE IN  -- shared coherence axis
# coherence from EDGES; >=2 regulators/targets ONLY (single-link genes trivially 1.0).
# Both panels SHOW the sign-flip null distribution (curve) instead of stating it in text.
from scipy import sparse as _sp
EDG=pd.concat([pd.read_csv(f"{FL}/flux_edges_{t}_CAD_aragam2022.tsv",sep="\t").assign(t=t) for t in TIS],ignore_index=True)
# IN: disease genes by their regulators ; OUT: broadcasters by their disease targets
inED=EDG[EDG.target_sig==True].assign(key=lambda d:d.t+"|"+d.target.map(base))
grpIN=inED.groupby("key")["flux"].apply(lambda v:np.array(v,float)); grpIN=grpIN[grpIN.apply(len)>=2]
grpOUT=EDG[EDG.target_sig==True].groupby([EDG.t,EDG.regulator.map(base)])["flux"].apply(lambda v:np.array(v,float))
grpOUT=grpOUT[grpOUT.apply(len)>10]                                # broadcasters (>10 disease targets)
bins=np.linspace(0,1,26); NPERM=200000; bsz=4000; NHIST=400

def coh_null(arrs,seed):
    """observed coherence + sign-flip null: count-null array (NPERM draws) + pooled null coherence sample."""
    obs_c=np.array([abs(v.sum())/np.abs(v).sum() for v in arrs]); ng=len(arrs)
    m=np.concatenate([np.abs(v) for v in arrs]); gi=np.concatenate([np.full(len(v),i) for i,v in enumerate(arrs)])
    ne=len(m); asum=np.array([np.abs(v).sum() for v in arrs])
    M=_sp.csr_matrix((np.ones(ne),(gi,np.arange(ne))),shape=(ng,ne))
    rng=np.random.default_rng(seed); cnt=[]; pool=None
    for s in range(0,NPERM,bsz):
        nb=min(bsz,NPERM-s); sg=rng.choice(np.array([-1.,1.]),size=(ne,nb))
        ck=np.abs(M@(sg*m[:,None]))/asum[:,None]                   # ng x nb coherence draws
        cnt.append((ck>THR).sum(0))
        if pool is None: pool=ck[:,:NHIST].ravel()                # sample for the null curve
    cnt=np.concatenate(cnt); z=(int((obs_c>THR).sum())-cnt.mean())/cnt.std()
    exc=int((cnt>=int((obs_c>THR).sum())).sum())
    pstr=f"p = {(exc+1)/(NPERM+1):.0e}" if exc else f"p < {1/(NPERM+1):.0e}"
    return obs_c,cnt,pool,z,pstr,ng

def cached_coh_null(key,arrs,seed):
    """Compute the sign-flip null once and cache it; reuse on later renders.
    Set env RECOMPUTE_NULL=1 (or delete the .npz) to force a fresh permutation."""
    cache=f"{FL}/perm_null_{key}_CAD_aragam2022.npz"
    if os.path.exists(cache) and os.environ.get("RECOMPUTE_NULL","0")!="1":
        d=np.load(cache,allow_pickle=True)
        print(f"  loaded cached null  {cache}")
        return d["obs_c"],d["cnt"],float(d["z"]),str(d["pstr"]),int(d["ng"])
    obs_c,cnt,_,z,pstr,ng=coh_null(arrs,seed)
    np.savez(cache,obs_c=obs_c,cnt=cnt,z=z,pstr=pstr,ng=ng,NPERM=NPERM,THR=THR)
    print(f"  saved null -> {cache}")
    return obs_c,cnt,z,pstr,ng

def panel(ax,arrs,seed,col,xlab,title,sig,show_n=True,key=None,legloc="upper right"):
    # permutation test: null distribution of the % of genes that are coherent (>0.7),
    # with the OBSERVED % drawn as a line -> reader sees where real data lands vs null.
    obs_c,cnt,z,pstr,ng=cached_coh_null(key,arrs,seed)
    obs=int((obs_c>THR).sum()); pn=100.0*cnt/ng; po=100.0*obs/ng; mu=pn.mean()
    ax.hist(pn,bins=45,color="#cdcdcd",edgecolor="white",lw=.3,density=True,label="sign-flip null",zorder=2)
    ax.axvline(mu,color="#777",ls="--",lw=1.1,zorder=3)
    ax.axvline(po,color=col,lw=2.8,zorder=4,label=f"observed  ({po:.0f}%)")
    lo=min(pn.min(),po); hi=max(pn.max(),po); pad=(hi-lo)*0.16+0.6
    ax.set_xlim(lo-pad,hi+pad); ax.set_ylim(0,ax.get_ylim()[1]*1.25); top=ax.get_ylim()[1]
    ax.set_xlabel(xlab); ax.set_ylabel("null density"); ax.set_yticks([])
    ax.legend(frameon=False,fontsize=7.6,loc=legloc)
    if show_n:
        ax.text(0.98,0.02,f"{ng} genes · {NPERM//1000}k sign-flip perm",transform=ax.transAxes,
                ha="right",va="bottom",fontsize=6.5,color="#999")
    ax.set_title(title,loc="left",fontweight="bold",fontsize=9.7)
    return z

# (c) KEY DRIVERS  &  (d) CORE MASTER REGULATORS -- own signal vs broadcasting breadth.
# Both panels share the SAME own-signal axis |z_TWAS| (a z-score): key drivers sit entirely
# BELOW the FDR line (non-significant), master regulators entirely ABOVE it (significant) --
# the differing regimes are the point, so a common y-scale makes them directly comparable.
# The x-axis (out-flux) is the total flux a regulator broadcasts to its disease targets.
N["bcast2"]=N.n_disease_targets>10
kd =N[(~N.sig)&N.bcast2].copy()                          # key drivers: non-significant broadcasters
cmr=N[N.sig & N.bcast2].copy()                           # core master regulators: significant broadcasters
kdl=kd[(kd.n_disease_targets>20)&(kd.out_flux>5)]        # key drivers worth labelling
zmax=ZMAX                                                # shared own-signal (|z_TWAS|) ceiling (same as panel a)
ZLAB="|z$_{\\mathrm{TWAS}}$|  (own signal, z-score)"
XLAB="out-flux  (total flux broadcast to disease genes)"

# (c) key drivers -- own axis to 5 (all key drivers sit below the FDR line), label ALL of them
C.scatter(kd.out_flux,kd.az,s=22,c="#0d3a66",alpha=0.9,edgecolors="white",linewidths=0.3,zorder=3)
C.set_ylim(0,3.5); C.set_yticks([0,1,2,3])
C.axhline(thr,ls="--",lw=1,color="#888"); C.text(C.get_xlim()[1],thr+0.05,"FDR 0.05",ha="right",va="bottom",fontsize=7.5,color="#888")
tC=[C.text(r.out_flux,r.az,f"{r.sym} ({r.tissue}, {int(r.n_disease_targets)})",fontsize=4.8,color="#0d3a66",fontstyle="italic") for _,r in kd.iterrows()]
adjust_text(tC,ax=C,arrowprops=dict(arrowstyle="-",color="#bbb",lw=0.35),expand=(1.05,1.25),force_text=(0.35,0.5))
C.set_xlabel(XLAB); C.set_ylabel(ZLAB)
C.set_title(f"(c) Key drivers — own signal vs broadcast  (n = {len(kd)})",loc="left",fontweight="bold",fontsize=9.7)

# (d) core master regulators -- label all
D.scatter(cmr.out_flux,cmr.az,s=42,c="#c1121f",alpha=0.9,edgecolors="white",linewidths=0.3,zorder=3)
D.set_ylim(0,zmax)
D.axhline(thr,ls="--",lw=1,color="#888"); D.text(D.get_xlim()[1],thr+0.15,"FDR 0.05",ha="right",va="bottom",fontsize=7.5,color="#888")
tD=[D.text(r.out_flux,r.az,f"{r.sym} ({r.tissue}, {int(r.n_disease_targets)})",fontsize=6.2,color="#7a0e15",fontstyle="italic") for _,r in cmr.iterrows()]
adjust_text(tD,ax=D,arrowprops=dict(arrowstyle="-",color="#bbb",lw=0.4),expand=(1.05,1.3),force_text=(0.4,0.6))
D.set_xlabel(XLAB); D.set_ylabel(ZLAB)
D.set_title(f"(d) Core master regulators — own signal vs broadcast  (n = {len(cmr)})",loc="left",fontweight="bold",fontsize=9.7)

# (e) COORDINATED OUT  &  (f) CONVERGENCE IN -- permutation-null panels
panel(E,list(grpOUT),1,"#e08214",
      "broadcasters coherent  (out-coherence > 0.7,  %)",
      "(e) Coordinated OUT — within null (incoherent)",False,show_n=False,key="out")
zF=panel(F,list(grpIN),0,"#2166ac",
      "disease genes coherent  (in-coherence > 0.7,  %)",
      "(f) Convergence IN — far beyond null",True,show_n=False,key="in",legloc="upper left")

for a in [A,B,E,F]: a.spines[["top","right"]].set_visible(False); a.grid(axis="y",color="0.93",lw=.5)
for a in [C,D]: a.spines[["top","right"]].set_visible(False); a.grid(axis="both",color="0.93",lw=.5)
fig.tight_layout()
for ext in ("pdf","png"): fig.savefig(f"{OUTDIR}/fig_main_readouts.{ext}",bbox_inches="tight")
print(f"panels e/f: out-broadcasters n={len(grpOUT)} ; in-disease-genes(>=2 reg) n={len(grpIN)} ; convergence z={zF:.0f}")
print(f"wrote {OUTDIR}/fig_main_readouts.pdf")
