#!/usr/bin/env python3
"""Data-leakage control: compare the full-cohort run with the 80%-train re-run.

Using the same individuals to reconstruct the network and to train the predictors
could, in principle, let full-cohort feature selection inflate downstream
associations. As a nested control we re-ran the entire procedure (network
reconstruction, channel-aware predictor training, and the TWAS LD/variance term)
on a random 80% training subset per tissue, so held-out samples never enter any
step. This script compares the two runs on the population-level quantities our
conclusions rest on, and writes a supplementary figure.

Quantities compared
    - significant genes per tissue and pooled (FDR < 0.05)
    - trans-using fraction among predictable two-channel genes (prediction step)
    - z_TWAS concordance for genes significant in both runs (Pearson r, sign
      agreement), gene-level on each gene's most-significant tissue
    - overlap of the significant gene sets

Outputs
    fig_leakage_traingrn.{png,pdf}   3-panel figure, in the shared figures directory
    leakage_per_tissue.tsv           per-tissue counts and trans-using fractions, in --out-dir
    leakage_summary.tsv              pooled headline numbers, in --out-dir

Input runs resolve under RESULTS (default results_cv/, override with FLUX_RESULTS); the figure
goes to results/figures/ by default (override with FLUX_FIGURES).
"""
from __future__ import annotations
import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
})

# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS = Path(os.environ.get("FLUX_RESULTS", Path(__file__).resolve().parents[2] / "results_cv"))
# The figure goes to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
FIGDIR = Path(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))
TISSUES = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
FDR = 0.05
PRED_R2 = 0.01
TRANS_CONFIGS = ("cis_trans", "trans_only")
FULL_BLUE = "#2C6FB5"
TRAIN_AMBER = "#E1812C"


def load_association(model_dir: Path, trait: str) -> pd.DataFrame:
    files = [f for f in glob.glob(str(model_dir / f"association_*_{trait}.tsv"))
             if "summary" not in os.path.basename(f)]
    return pd.concat([pd.read_csv(f, sep="\t") for f in files], ignore_index=True)


def load_selected(model_dir: Path) -> pd.DataFrame:
    files = glob.glob(str(model_dir / "*.horseshoe_selected.tsv"))
    return pd.concat(
        [pd.read_csv(f, sep="\t").assign(tissue=os.path.basename(f).split(".")[0])
         for f in files],
        ignore_index=True,
    )


def trans_using_fraction(selected: pd.DataFrame) -> pd.Series:
    """Per-tissue trans-using fraction among predictable two-channel genes.

    A gene counts only if it is predictable (LOO-R2 >= PRED_R2) and has both a cis
    and a trans channel available, so the fraction measures how often the model
    selection prefers trans when it genuinely has the choice.
    """
    pred = selected[selected.best_loo_r2 >= PRED_R2]
    two = pred[(pred.n_cis > 0) & (pred.n_trans > 0)].copy()
    two["trans"] = two.best_config.isin(TRANS_CONFIGS)
    return two.groupby("tissue").trans.mean()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full-dir", default=str(RESULTS / "horseshoe_alltargets"))
    ap.add_argument("--train-dir", default=str(RESULTS / "horseshoe_alltargets_traingrn"))
    ap.add_argument("--trait", default="CAD_aragam2022")
    ap.add_argument("--out-dir", default=str(RESULTS / "leakage_check"))
    ap.add_argument("--fig-dir", default=str(FIGDIR),
                    help="Directory for the supplementary figure (default: the shared "
                         "figures directory).")
    args = ap.parse_args()

    full_dir, train_dir = Path(args.full_dir), Path(args.train_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    assoc_full = load_association(full_dir / "association", args.trait)
    assoc_train = load_association(train_dir / "association", args.trait)
    sel_full = load_selected(full_dir / "models")
    sel_train = load_selected(train_dir / "models")

    sig_full = assoc_full[assoc_full.p_adj < FDR]
    sig_train = assoc_train[assoc_train.p_adj < FDR]

    # --- per-tissue significant counts and trans-using fractions ---
    cnt_full = sig_full.groupby("tissue").gene_id.size()
    cnt_train = sig_train.groupby("tissue").gene_id.size()
    tu_full = trans_using_fraction(sel_full)
    tu_train = trans_using_fraction(sel_train)

    per_tissue = pd.DataFrame({
        "tissue": TISSUES,
        "sig_full": [int(cnt_full.get(t, 0)) for t in TISSUES],
        "sig_train": [int(cnt_train.get(t, 0)) for t in TISSUES],
        "trans_using_full": [tu_full.get(t, np.nan) for t in TISSUES],
        "trans_using_train": [tu_train.get(t, np.nan) for t in TISSUES],
    })
    per_tissue.to_csv(out_dir / "leakage_per_tissue.tsv", sep="\t", index=False)

    # --- pooled headline numbers ---
    genes_full = set(sig_full.gene_id)
    genes_train = set(sig_train.gene_id)
    shared = genes_full & genes_train

    pred_full = sel_full[sel_full.best_loo_r2 >= PRED_R2]
    pred_train = sel_train[sel_train.best_loo_r2 >= PRED_R2]
    two_full = pred_full[(pred_full.n_cis > 0) & (pred_full.n_trans > 0)]
    two_train = pred_train[(pred_train.n_cis > 0) & (pred_train.n_trans > 0)]
    tu_pooled_full = two_full.best_config.isin(TRANS_CONFIGS).mean()
    tu_pooled_train = two_train.best_config.isin(TRANS_CONFIGS).mean()

    # gene-level concordance: genes significant in both, most-significant tissue
    def best_tissue_z(assoc: pd.DataFrame) -> pd.DataFrame:
        sub = assoc[assoc.gene_id.isin(shared)]
        keep = sub.loc[sub.groupby("gene_id").p_adj.idxmin()]
        return keep[["gene_id", "z_twas"]]

    merged = best_tissue_z(assoc_full).merge(
        best_tissue_z(assoc_train), on="gene_id", suffixes=("_full", "_train"))
    r = float(np.corrcoef(merged.z_twas_full, merged.z_twas_train)[0, 1])
    sign_agree = float((np.sign(merged.z_twas_full)
                        == np.sign(merged.z_twas_train)).mean())

    summary = pd.DataFrame([
        ("sig_genes_full", len(genes_full)),
        ("sig_genes_train", len(genes_train)),
        ("sig_assoc_full", len(sig_full)),
        ("sig_assoc_train", len(sig_train)),
        ("shared_genes", len(shared)),
        ("overlap_shared_over_full", len(shared) / len(genes_full)),
        ("trans_using_full", tu_pooled_full),
        ("trans_using_train", tu_pooled_train),
        ("concordance_pearson_r", r),
        ("concordance_sign_agreement", sign_agree),
        ("concordance_n_genes", len(merged)),
    ], columns=["metric", "value"])
    summary.to_csv(out_dir / "leakage_summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))
    x = np.arange(len(TISSUES))
    w = 0.4

    # (a) significant genes per tissue
    ax = axes[0]
    ax.bar(x - w / 2, per_tissue.sig_full, w, color=FULL_BLUE, label="Full cohort")
    ax.bar(x + w / 2, per_tissue.sig_train, w, color=TRAIN_AMBER, label="80% train")
    ax.set_xticks(x)
    ax.set_xticklabels(TISSUES, rotation=45, ha="right")
    ax.set_ylabel("Significant genes (FDR $<$ 0.05)")
    ax.set_title("(a) Associations per tissue", fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=8)

    # (b) trans-using fraction per tissue
    ax = axes[1]
    ax.bar(x - w / 2, per_tissue.trans_using_full * 100, w, color=FULL_BLUE)
    ax.bar(x + w / 2, per_tissue.trans_using_train * 100, w, color=TRAIN_AMBER)
    ax.set_xticks(x)
    ax.set_xticklabels(TISSUES, rotation=45, ha="right")
    ax.set_ylabel("Trans-using among two-channel\npredictable genes (%)")
    ax.set_title("(b) Trans usage per tissue", fontsize=9, loc="left")

    # (c) gene-level z concordance
    ax = axes[2]
    lim = np.ceil(np.abs(merged[["z_twas_full", "z_twas_train"]].values).max() / 5) * 5
    ax.axhline(0, color="0.8", lw=0.6)
    ax.axvline(0, color="0.8", lw=0.6)
    ax.plot([-lim, lim], [-lim, lim], ls="--", color="0.5", lw=0.8)
    ax.scatter(merged.z_twas_full, merged.z_twas_train, s=8,
               color=FULL_BLUE, alpha=0.45, edgecolors="none")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("$z_{\\mathrm{TWAS}}$ (full cohort)")
    ax.set_ylabel("$z_{\\mathrm{TWAS}}$ (80% train)")
    ax.set_title("(c) Effect concordance", fontsize=9, loc="left")
    ax.text(0.04, 0.96,
            f"$r = {r:.2f}$\n{sign_agree*100:.0f}% sign agree\n$n = {len(merged)}$",
            transform=ax.transAxes, va="top", ha="left", fontsize=8)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_leakage_traingrn.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"[leakage] wrote tables to {out_dir}/ and figure to {fig_dir}/")


if __name__ == "__main__":
    main()
