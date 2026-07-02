#!/usr/bin/env python3
"""
Shared data layer + style for the reconstructed GRN-TWAS figures.

Centralises three things every figure needs so the panels stay mutually
consistent:

  * STYLE / palette   -- apply_style(), TISSUES, TISSUE_COLORS, the cis /
    cis+trans / trans-only category palette, and resolve_symbol (HGNC labels).
  * GRN best-of model  -- load_best_of(): per (tissue, gene_id), the held-out
    R^2 of each channel taken as the BEST across all Bayesian methods present
    (bayes_ridge / bvsr / horseshoe / bslmm, pooled over one or more predictor
    dirs), then each predictable gene assigned to ONE non-overlapping category
    (cis_only / cis_trans / trans_only) via the same rule as
    plot_grn_vs_classical.py.
  * feature counts     -- load_feature_counts(): per gene, #cis (Stage 4
    LD-pruned), #trans (Stage 7 GRN-walk SNPs) and #GRN parents (Stage 7 hop-1
    source genes).

Data-agnostic: all paths are passed in. Gene symbols come from gene_labels.py
(GENCODE v19); pass the annot path or set $GENE_ANNOT.
"""
from __future__ import annotations

import glob
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# canonical STARNET tissue order + a stable per-tissue colour
TISSUES = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
TISSUE_COLORS = {
    "AOR": "#e41a1c", "Blood": "#377eb8", "LIV": "#4daf4a", "MAM": "#984ea3",
    "SF": "#ff7f00", "SKLM": "#a65628", "VAF": "#f781bf",
}

# three mutually-exclusive predictor categories (hop depth not distinguished)
CATS = ["cis_only", "cis_trans", "trans_only"]
CAT_COLORS = {"cis_only": "#2166ac", "cis_trans": "#1b7837", "trans_only": "#d73027"}
CAT_LABELS = {"cis_only": "Cis", "cis_trans": "Cis + trans", "trans_only": "Trans only"}

# config names in the Stage 9 scores -> the three channels
CIS, TRANS, JOINT = "cis_only", "trans_only", "cis_trans_resid"


def apply_style(base: int = 9) -> None:
    """Publication style: Arial, no top/right spines, 42-type fonts for vector edit."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base, "axes.titlesize": base + 1, "axes.labelsize": base,
        "xtick.labelsize": base - 1, "ytick.labelsize": base - 1,
        "legend.fontsize": base - 0.5, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


# --- gene symbols -----------------------------------------------------------

def _gene_labels():
    import importlib.util
    here = os.path.join(os.path.dirname(os.path.realpath(__file__)), "gene_labels.py")
    spec = importlib.util.spec_from_file_location("gene_labels", here)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def symbol_map(gene_annot: str | None = None) -> dict:
    gl = _gene_labels()
    return gl.load_label_map(gl.default_annot(gene_annot))


# --- GRN best-of-method predictability --------------------------------------

def _read_scores(pred_dirs) -> pd.DataFrame:
    """Concatenate every <tissue>.bayesian_scores.tsv under the given dir(s),
    tagging each row with its tissue (parsed from the file basename)."""
    parts = []
    for d in pred_dirs:
        for p in sorted(glob.glob(os.path.join(d, "*.bayesian_scores.tsv"))):
            t = os.path.basename(p).split(".")[0]
            df = pd.read_csv(p, sep="\t", na_values=["NA"])
            df["tissue"] = t
            parts.append(df)
    if not parts:
        raise SystemExit(f"no *.bayesian_scores.tsv under {pred_dirs}")
    return pd.concat(parts, ignore_index=True)


def classify(out: pd.DataFrame, thr: float, delta: float) -> pd.DataFrame:
    """Non-overlapping channel per predictable gene (verbatim from
    plot_grn_vs_classical.py.classify, GRN side only)."""
    cis = out["r2_cis"].fillna(-np.inf)
    trans = out["r2_trans"].fillna(-np.inf)
    joint = out["r2_joint"].fillna(-np.inf)
    grn_pred = out["r2_grn"] >= thr
    cat = np.where(~grn_pred, "none", "")
    m_trans = grn_pred & (cis < thr) & (trans >= thr)
    m_cistrans = grn_pred & ~m_trans & (
        ((cis < thr) & (joint >= thr)) | ((cis >= thr) & ((joint - cis) > delta)))
    m_cis = grn_pred & ~m_trans & ~m_cistrans
    cat = np.where(m_trans, "trans_only", cat)
    cat = np.where(m_cistrans, "cis_trans", cat)
    cat = np.where(m_cis, "cis_only", cat)
    out["grn_cat"] = cat
    return out


def load_best_of(pred_dirs, thr: float = 0.01, delta: float = 0.01,
                 methods=None) -> pd.DataFrame:
    """Per (tissue, gene_id): the GRN-TWAS predictor = the single (method, config)
    with the highest held-out R^2 -- the per-gene ARGMAX over the Bayesian methods
    x channels (exactly what Stage 10 does). Because every method's own cis_only is
    a candidate, a `cis_trans` winner NECESSARILY beat that method's cis-only, so
    the trans component genuinely ADDS over cis WITHIN the chosen model (a true
    cis+trans gene -- never a cis gene whose trans makes it worse). The winning
    model's three channel R^2 are returned as r2_cis / r2_trans / r2_joint (so the
    within-model cis baseline and trans gain are self-consistent).
    Columns: r2_cis, r2_trans, r2_joint, r2_grn (=best), win_method, win_config,
    grn_cat in {cis_only, cis_trans, trans_only, none}."""
    bay = _read_scores(pred_dirs)
    bay = bay[bay["config"].isin([CIS, TRANS, JOINT])].dropna(subset=["test_r2"])
    if methods is not None:                       # restrict to a single method (e.g. bslmm)
        bay = bay[bay["method"].isin(methods)]
    win_idx = bay.groupby(["tissue", "gene_id"])["test_r2"].idxmax()
    win = (bay.loc[win_idx, ["tissue", "gene_id", "method", "config", "test_r2"]]
              .rename(columns={"method": "win_method", "config": "win_config",
                               "test_r2": "r2_grn"}))
    # the winning method's own cis_only / trans_only / cis_trans_resid R^2
    wm = bay.merge(win[["tissue", "gene_id", "win_method"]], on=["tissue", "gene_id"])
    wm = wm[wm["method"] == wm["win_method"]]
    chan = (wm.pivot_table(index=["tissue", "gene_id"], columns="config",
                           values="test_r2", aggfunc="first").reset_index())
    for c in (CIS, TRANS, JOINT):
        if c not in chan.columns:
            chan[c] = np.nan
    chan = chan.rename(columns={CIS: "r2_cis", TRANS: "r2_trans", JOINT: "r2_joint"})
    out = win.merge(chan, on=["tissue", "gene_id"])
    out["grn_cat"] = out["win_config"].map(
        {CIS: "cis_only", JOINT: "cis_trans", TRANS: "trans_only"})
    out.loc[out["r2_grn"] < thr, "grn_cat"] = "none"
    return out


# --- feature counts (Stage 4 cis, Stage 7 trans + GRN parents) --------------

def load_feature_counts(ld_prune_dir: str, trans_dir: str, tissue: str) -> pd.DataFrame:
    """Per gene_id for one tissue: n_cis (LD-pruned cis SNPs), n_trans (unique
    GRN-walk trans SNPs) and n_parents (distinct hop-1 GRN source genes)."""
    cisp = Path(ld_prune_dir) / f"{tissue}.cis_snps.pruned.tsv.gz"
    trp = Path(trans_dir) / f"{tissue}.trans_features.tsv.gz"
    n_cis = (pd.read_csv(cisp, sep="\t", usecols=["gene_id", "rs_id"])
               .groupby("gene_id")["rs_id"].nunique().rename("n_cis"))
    tr = pd.read_csv(trp, sep="\t", usecols=["gene_id", "hop", "source_gene_id", "rs_id"])
    n_trans = tr.groupby("gene_id")["rs_id"].nunique().rename("n_trans")
    n_par = (tr[tr["hop"] == 1].groupby("gene_id")["source_gene_id"]
               .nunique().rename("n_parents"))
    out = pd.concat([n_cis, n_trans, n_par], axis=1)
    out["n_cis"] = out["n_cis"].fillna(0).astype(int)
    out["n_trans"] = out["n_trans"].fillna(0).astype(int)
    out["n_parents"] = out["n_parents"].fillna(0).astype(int)
    out["tissue"] = tissue
    return out.reset_index()
