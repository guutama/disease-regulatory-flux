#!/usr/bin/env python3
"""Cardiovascular functional test and universe builder (results_findr_py).

Builds the curated cardiovascular gene universe from CV-matching terms across five downloaded
Enrichr libraries (MGI mouse-knockout phenotypes, GO Biological Process, Reactome, KEGG, Jensen
DISEASES; data/genesets/), writes the universe and the selected terms out as data, and tests
whether the disease-gene predictor classes (cis-only, cis+trans, trans-only) and the convergent
set are enriched for cardiovascular function over the tested-gene background (one-sided Fisher).
MGI mouse-knockout phenotypes are a readout functionally orthogonal to human genetics.

All inputs and outputs are under results_findr_py/, so the universe behind the functional figure is
a saved, reproducible artifact from the same pipeline as the rest of the paper.

Outputs (results_findr_py/functional/):
  cardiovascular_terms.tsv      one row per selected cardiovascular term: library, term, n_genes
  cardiovascular_universe.tsv   one row per cardiovascular gene: gene, in_tested, n_cv_terms,
                                mgi_ko_cv, example_terms (the union, before/after intersecting tested)
  cardiovascular_enrichment.tsv one row per gene class x axis: enrichment odds ratio, 95% CI,
                                one-sided Fisher p-value and overlap counts (the figure's numbers)
"""
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "utils"))
import gene_labels as gl

RF = ROOT / "results_findr_py"
AT = str(RF)
TIS = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
TRAIT = "CAD_aragam2022"
GS = str(ROOT / "data" / "genesets")
OUTD = str(RF / "functional")
MAPF = os.environ.get(
    "GENE_ANNOT",
    "/cluster/projects/nn1015k/datasets/references/gencode/GRCh37/gencode.v19.genes.tsv")
base = lambda g: str(g).split(".")[0]
MAP = gl.load_label_map(MAPF)
sym = lambda g: gl.resolve(g, MAP)
CV = re.compile(
    r"cardiovascular|cardiac|heart|myocard|cardiomyo|vascul|vasocon|vasodil|"
    r"angiogen|endotheli|blood vessel|artery|arterial|aorta|aortic|atheroscler|"
    r"smooth muscle|blood pressure|hypertens|coagul|platelet|thrombo|hemostasis|fibrin|"
    r"lipid|lipoprotein|cholesterol|choleste|triglycer|sterol|foam cell|"
    r"coronary|ischemi|infarct|elastic fiber|extracellular matrix", re.I)

# Five downloaded Enrichr libraries; readable library label -> filename.
GMT = [("MGI", "MGI_Mammalian_Phenotype_Level_4_2021.gmt"),
       ("GO-BP", "GO_Biological_Process_2023.gmt"),
       ("Reactome", "Reactome_2022.gmt"),
       ("KEGG", "KEGG_2021_Human.gmt"),
       ("Jensen", "Jensen_DISEASES.gmt")]


def load_gmt(path):
    d = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 3:
            d[f[0]] = set(g.upper() for g in f[2:] if g)
    return d


# CV-matching terms, keyed by (library, term) so identical term names in different libraries
# do not collide (matching the figure's library-prefixed universe).
CVTERMS = {}
for lib, fn in GMT:
    for t, gs in load_gmt(f"{GS}/{fn}").items():
        if CV.search(t):
            CVTERMS[(lib, t)] = gs
print(f"CV terms collected: {len(CVTERMS)} across {len(GMT)} libraries", flush=True)

# gene classes from the association tables (selected-config predictor per gene)
A = pd.concat([pd.read_csv(f"{AT}/association/association_{t}_{TRAIT}.tsv", sep="\t") for t in TIS],
              ignore_index=True)
if "win_config" not in A.columns and "config" in A.columns:
    A["win_config"] = A["config"]
A["b"] = A.gene_id.map(base)
A = A.sort_values("p_adj")
best = A.drop_duplicates("b", keep="first")
tested = set(A.b.unique())
sig = best[best.p_adj < 0.05]
S = pd.read_csv(f"{AT}/functional/functional_relevance_summary.tsv", sep="\t")
conv = set(S.loc[S.convergent, "tbase"])
nonconv = set(S.loc[~S.convergent, "tbase"])
cisonly = set(sig.loc[sig.win_config == "cis_only", "b"])
# trans-using: significant genes whose selected model uses the trans channel (cis+trans or trans-only)
transusing = set(sig.loc[sig.win_config.isin(["cis_trans", "trans_only"]), "b"])
conv_TU = conv & transusing
tested_sym = set(sym(g) for g in tested)

# non-significant contributors: genes that feed trans signal into disease genes through the network
# but are not themselves disease-significant (the regulators of the flux map), pooled to unique genes.
NODES = str(RF / "flux" / "nodes" / "flux_nodes_{t}_CAD_aragam2022.tsv")
_tru = lambda x: str(x).strip().lower() in ("true", "1")
noncontrib = set()
for t in TIS:
    nd = pd.read_csv(NODES.format(t=t), sep="\t")
    nd["b"] = nd.gene_id.map(base)
    noncontrib |= set(nd.loc[~nd.disease_sig.map(_tru), "b"])

# ---- write the cardiovascular gene sets (terms) and the universe as data ----
os.makedirs(OUTD, exist_ok=True)
terms_rows = [{"library": lib, "term": t, "n_genes": len(gs)} for (lib, t), gs in CVTERMS.items()]
pd.DataFrame(terms_rows).sort_values(["library", "term"]).to_csv(
    f"{OUTD}/cardiovascular_terms.tsv", sep="\t", index=False)

mgi_genes = set().union(*[gs for (lib, t), gs in CVTERMS.items() if lib == "MGI"])
allcv = set().union(*CVTERMS.values())
uni_rows = []
for g in sorted(allcv):
    hit = [f"{lib}:{t}" for (lib, t), gs in CVTERMS.items() if g in gs]
    uni_rows.append({"gene": g, "in_tested": g in tested_sym, "n_cv_terms": len(hit),
                     "mgi_ko_cv": g in mgi_genes, "example_terms": "; ".join(hit[:3])})
pd.DataFrame(uni_rows).to_csv(f"{OUTD}/cardiovascular_universe.tsv", sep="\t", index=False)

CVgenes = allcv & tested_sym
MGIcv = mgi_genes & tested_sym
print(f"wrote {OUTD}/cardiovascular_terms.tsv ({len(terms_rows)} terms) and "
      f"{OUTD}/cardiovascular_universe.tsv ({len(uni_rows)} genes)")
print(f"CV gene universe (in tested set): {len(CVgenes)} of {len(tested_sym)} tested "
      f"({100*len(CVgenes)/len(tested_sym):.0f}%); MGI-CV {len(MGIcv)} "
      f"({100*len(MGIcv)/len(tested_sym):.0f}%)\n")


TESTED_LIST = sorted(tested_sym)
RNG = np.random.default_rng(0)
NPERM = 5000                                              # random size-matched null draws per class


def in_cv_test(geneset, universe, label):
    qs = set(sym(g) for g in geneset) & tested_sym
    n = len(qs)
    a = len(qs & universe)
    rest = tested_sym - qs
    c = len(rest & universe)
    OR, p = stats.fisher_exact([[a, n - a], [c, len(rest) - c]], alternative="greater")
    aa, bb, cc, dd = (a + 0.5, n - a + 0.5, c + 0.5, len(rest) - c + 0.5)
    orr = (aa * dd) / (bb * cc)
    se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    # random size-matched null: OR of NPERM random n-gene sets drawn from the tested background
    cvmask = np.fromiter((g in universe for g in TESTED_LIST), bool, len(TESTED_LIST))
    N = len(TESTED_LIST)
    Ncv = int(cvmask.sum())
    nulls = np.empty(NPERM)
    for b in range(NPERM):
        k = int(cvmask[RNG.choice(N, size=n, replace=False)].sum())
        na, nb, nc, nd = (k + 0.5, n - k + 0.5, Ncv - k + 0.5, (N - n) - (Ncv - k) + 0.5)
        nulls[b] = (na * nd) / (nb * nc)
    perm_p = (int(np.sum(nulls >= orr)) + 1) / (NPERM + 1)
    print(f"  {label:28s} n={n:4d}  in-CV {100*a/max(n,1):4.1f}%  OR={OR:.2f}  p={p:.2g}  "
          f"null_OR={nulls.mean():.2f}  perm_p={perm_p:.1e}")
    return dict(n_genes=n, n_in_universe=a, odds_ratio=round(orr, 3),
                ci_low=round(float(np.exp(np.log(orr) - 1.96 * se)), 3),
                ci_high=round(float(np.exp(np.log(orr) + 1.96 * se)), 3),
                p_value=p, significant=bool(p < 0.05),
                null_or_mean=round(float(nulls.mean()), 3),
                null_or_ci_low=round(float(np.percentile(nulls, 2.5)), 3),
                null_or_ci_high=round(float(np.percentile(nulls, 97.5)), 3),
                perm_p=perm_p)


# Four classes: three of disease genes (by how their association is built) plus the non-significant
# contributors (the regulators the network implicates). Convergence is a single class.
CLASSES = [("cis_only", cisonly), ("trans_using", transusing), ("convergent", conv),
           ("non-significant contributor", noncontrib)]
erows = []
print("=== enrichment in cardiovascular gene universe (vs tested background) ===")
for lab, gs in CLASSES:
    r = in_cv_test(gs, CVgenes, lab)
    erows.append({"gene_class": lab, "axis": "cardiovascular_universe", **r})
print(f"\n=== MGI mouse-KO cardiovascular-phenotype genes only ({len(MGIcv)} in tested) ===")
for lab, gs in CLASSES:
    r = in_cv_test(gs, MGIcv, lab)
    erows.append({"gene_class": lab, "axis": "mgi_knockout_cardiovascular", **r})
pd.DataFrame(erows).to_csv(f"{OUTD}/cardiovascular_enrichment.tsv", sep="\t", index=False)
print(f"\nwrote {OUTD}/cardiovascular_enrichment.tsv ({len(erows)} class x axis rows)")


def top_terms(geneset, topn=8, kmin=3, Kmin=5, Kmax=1500):
    """Cardiovascular terms most over-represented in a class, by fold enrichment (observed/expected)."""
    qs = set(sym(g) for g in geneset) & tested_sym
    n = len(qs)
    N = len(tested_sym)
    rows = []
    for (lib, term), gs in CVTERMS.items():
        gsT = gs & tested_sym
        K = len(gsT)
        if K < Kmin or K > Kmax:
            continue
        k = len(gsT & qs)
        if k < kmin:
            continue
        rows.append({"library": lib, "term": term, "genes_hit": k, "term_size": K,
                     "fold_enrichment": round(k / (n * K / N), 3)})
    return pd.DataFrame(rows).sort_values("fold_enrichment", ascending=False).head(topn)


trows = []
for lab, gs in CLASSES:
    for _, r in top_terms(gs).iterrows():
        trows.append({"gene_class": lab, **r.to_dict()})
pd.DataFrame(trows).to_csv(f"{OUTD}/cardiovascular_top_terms.tsv", sep="\t", index=False)
print(f"wrote {OUTD}/cardiovascular_top_terms.tsv ({len(trows)} rows)")
