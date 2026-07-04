#!/usr/bin/env python3
"""Covariate-structure quality-control figures for the STARNET expression matrices.

The STARNET expression matrices were covariate-adjusted upstream (an R pipeline whose steps
include a linear-model removal of age, sex and genotype PCs). This script tests that claim
directly: for each tissue it computes the leading principal components of the samples x genes
expression matrix and measures each PC's association (Pearson correlation) with sample
covariates from the cohort workbook. If the adjustment worked, the leading components should
carry no age or sex structure.

Figures (supplementary)
    fig_covariate_check.{pdf,png}         per-tissue heatmap of -log10 P (PC x covariate),
                                          cyan stars marking Bonferroni-significant cells
    fig_covariate_pc_scatter.{pdf,png}    PC1-PC2 scatter of the representative tissue
                                          (default SKLM), coloured by age, sex and BMI

Quality-control tables (in --out-dir)
    pc_covariate_association_<TISSUE>.csv   per-tissue PC x covariate p-values
    pc_covariate_association_all.csv        long-format combined table
    fig_pc_scatter_<TISSUE>.{pdf,png}       PC1-PC2 scatter for every tissue
    SUMMARY.md                              human-readable conclusion

Only aggregate PC-covariate statistics and per-sample principal-component scatters are written;
no expression values or covariate values are exported. The cohort workbook is restricted STARNET
data and must be supplied explicitly with --pheno-xls (there is no default path).

Figures go to results/figures/ by default (override with FLUX_FIGURES). Style: Arial, 300 DPI.
"""
from __future__ import annotations
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
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

# All figures go to one shared directory: results/figures/ by default (override with FLUX_FIGURES).
FIGDIR = Path(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))

# Covariates to test. Numeric ones get Pearson r vs each PC; categorical/binary
# ones get a point-biserial / ANOVA association.  Age + Sex are the ones the
# upstream `lm` step claims to have regressed out.
NUMERIC_COVS = {
    "Age": "Age",
    "BMI": "BMI(kg/m2)",
    "SBP": "SBP",
    "DBP": "DBP",
    "Syntax_score": "Syntax sc",   # angiographic CAD severity
    "CRP": "CRP(mg/l)",
}
BINARY_COVS = {
    "Sex": ("Sex", {"male": 1.0, "female": 0.0}),
    "Smoker": ("Smoker", {"yes": 1.0, "no": 0.0}),
    "LipidLowerer": ("LipidLowerer", {"yes": 1.0, "no": 0.0}),
}

# Tissues reported in the primary analysis (FC and MP are excluded).
PAPER_TISSUES = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]

# Display labels used on figure axes (the internal keys above are terse).
COV_LABELS = {
    "Age": "Age",
    "BMI": "BMI",
    "SBP": "SBP",
    "DBP": "DBP",
    "Syntax_score": "Syntax score",
    "CRP": "CRP",
    "Sex": "Sex",
    "Smoker": "Smoker",
    "LipidLowerer": "Lipid-lowering",
}


def load_phenotypes(xls_path: Path) -> pd.DataFrame:
    df = pd.read_excel(xls_path, sheet_name="cohort")
    df["starnet-ID"] = df["starnet-ID"].astype(str)
    out = pd.DataFrame(index=df["starnet-ID"].values)
    out.index.name = "sample"
    for name, col in NUMERIC_COVS.items():
        if col in df.columns:
            out[name] = pd.to_numeric(df[col], errors="coerce").values
    for name, (col, mapping) in BINARY_COVS.items():
        if col in df.columns:
            raw = df[col].astype(str).str.strip().str.lower()
            out[name] = raw.map(mapping).values
    return out


def load_expression(path: Path) -> pd.DataFrame:
    """Return samples x genes matrix (index = sample id string)."""
    df = pd.read_csv(path, sep="\t")
    meta_cols = [c for c in ("gene_symbol", "gene_id") if c in df.columns]
    expr = df.drop(columns=meta_cols)
    # rows=genes, cols=samples -> transpose to samples x genes
    mat = expr.T
    mat.index = mat.index.astype(str)
    mat.index.name = "sample"
    return mat.apply(pd.to_numeric, errors="coerce")


def compute_pcs(mat: pd.DataFrame, n_pcs: int) -> tuple[pd.DataFrame, np.ndarray]:
    """PCA on samples x genes. Standardise genes, drop zero-variance genes."""
    X = mat.values.astype(float)
    # drop genes with any NaN or zero variance
    good = ~np.any(np.isnan(X), axis=0)
    Xg = X[:, good]
    sd = Xg.std(axis=0)
    Xg = Xg[:, sd > 0]
    Xg = (Xg - Xg.mean(axis=0)) / Xg.std(axis=0)
    # SVD on centred/scaled matrix
    U, S, _ = np.linalg.svd(Xg, full_matrices=False)
    k = min(n_pcs, S.shape[0])
    pcs = U[:, :k] * S[:k]
    var_explained = (S[:k] ** 2) / (S ** 2).sum()
    cols = [f"PC{i+1}" for i in range(k)]
    return pd.DataFrame(pcs, index=mat.index, columns=cols), var_explained


def associate(pc: np.ndarray, cov: np.ndarray) -> tuple[float, float]:
    """Pearson r and two-sided p between one PC and one covariate (NaN-safe)."""
    m = ~(np.isnan(pc) | np.isnan(cov))
    if m.sum() < 10 or np.nanstd(cov[m]) == 0:
        return np.nan, np.nan
    r, p = stats.pearsonr(pc[m], cov[m])
    return r, p


def var_for_pc(long: pd.DataFrame, tissue: str, pc: str) -> float:
    """Variance fraction explained by one PC of one tissue (for axis labels)."""
    sub = long[(long["tissue"] == tissue) & (long["pc"] == pc)]
    return float(sub["var_explained"].iloc[0]) if len(sub) else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--expr-dir", default="results/stage2_test/expression")
    ap.add_argument("--pheno-xls", required=True,
                    help="STARNET cohort workbook (.xls) with sample covariates. "
                         "Restricted data; must be supplied explicitly.")
    ap.add_argument("--out-dir", default="results/qc_covariate_check")
    ap.add_argument("--fig-dir", default=str(FIGDIR),
                    help="Directory for the supplementary figures (default: the shared "
                         "figures directory).")
    ap.add_argument("--rep-tissue", default="SKLM",
                    help="Tissue whose PC1-PC2 scatter becomes fig_covariate_pc_scatter "
                         "(default: SKLM, the tissue with the only residual associations).")
    ap.add_argument("--n-pcs", type=int, default=15)
    ap.add_argument(
        "--tissues",
        default=",".join(PAPER_TISSUES),
        help="comma-separated tissues to include (default: the seven "
             "tissues used in the primary analysis)",
    )
    args = ap.parse_args()
    keep_tissues = [t for t in args.tissues.split(",") if t]

    expr_dir = Path(args.expr_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(args.fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    ph = load_phenotypes(Path(args.pheno_xls))
    cov_names = list(ph.columns)
    print(f"[cov-test] phenotypes: {ph.shape[0]} samples, covariates: {cov_names}")

    all_files = sorted(expr_dir.glob("*.expr.tsv.gz"))
    expr_files = [f for f in all_files if f.name.split(".")[0] in keep_tissues]
    print(f"[cov-test] {len(all_files)} expression files present; "
          f"analysing {len(expr_files)}: "
          f"{', '.join(f.name.split('.')[0] for f in expr_files)}")

    all_rows = []
    per_tissue_pc_cov: dict[str, pd.DataFrame] = {}
    pcs_by_tissue: dict[str, pd.DataFrame] = {}
    aligned_ph: dict[str, pd.DataFrame] = {}

    for f in expr_files:
        tissue = f.name.split(".")[0]
        mat = load_expression(f)
        pcs, var_exp = compute_pcs(mat, args.n_pcs)
        # align phenotypes to this tissue's samples
        ph_t = ph.reindex(pcs.index)
        n_match = ph_t["Age"].notna().sum() if "Age" in ph_t else 0
        print(f"  {tissue}: {mat.shape[0]} samples x {mat.shape[1]} genes; "
              f"PC1 var={var_exp[0]:.3f}; matched-pheno={n_match}")

        grid = pd.DataFrame(index=pcs.columns, columns=cov_names, dtype=float)
        for cov in cov_names:
            cv = ph_t[cov].values.astype(float)
            for pc_name in pcs.columns:
                r, p = associate(pcs[pc_name].values, cv)
                grid.loc[pc_name, cov] = p
                all_rows.append({
                    "tissue": tissue, "pc": pc_name, "covariate": cov,
                    "var_explained": var_exp[int(pc_name[2:]) - 1],
                    "pearson_r": r, "p_value": p,
                })
        per_tissue_pc_cov[tissue] = grid
        pcs_by_tissue[tissue] = pcs
        aligned_ph[tissue] = ph_t
        grid.to_csv(out_dir / f"pc_covariate_association_{tissue}.csv")

    long = pd.DataFrame(all_rows)
    long.to_csv(out_dir / "pc_covariate_association_all.csv", index=False)

    # Bonferroni threshold across all tests
    n_tests = long["p_value"].notna().sum()
    bonf = 0.05 / max(n_tests, 1)
    print(f"[cov-test] {n_tests} tests; Bonferroni alpha = {bonf:.2e}")

    # ---- Figure 1 (fig_covariate_check): heatmap grid of -log10 p, one panel per tissue ----
    tissues = list(per_tissue_pc_cov.keys())
    cov_labels = [COV_LABELS.get(c, c) for c in cov_names]
    ncol = 3
    nrow = int(np.ceil(len(tissues) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.5 * nrow),
                             squeeze=False, constrained_layout=True)
    vmax = 10
    for ax in axes.flat:
        ax.set_visible(False)
    im = None
    for i, tissue in enumerate(tissues):
        ax = axes[i // ncol][i % ncol]
        ax.set_visible(True)
        grid = per_tissue_pc_cov[tissue].astype(float)
        neglogp = -np.log10(grid.values.clip(min=1e-300))
        im = ax.imshow(neglogp, aspect="auto", cmap="magma",
                       vmin=0, vmax=vmax)
        ax.set_xticks(range(len(cov_labels)))
        # only label the x-axis on the bottom panel of each column
        bottom_row = (i // ncol) == (nrow - 1) or (i + ncol) >= len(tissues)
        if bottom_row:
            ax.set_xticklabels(cov_labels, rotation=90, fontsize=7)
        else:
            ax.set_xticklabels([])
        ax.set_yticks(range(grid.shape[0]))
        ax.set_yticklabels(grid.index, fontsize=6)
        ax.tick_params(length=0)
        panel = chr(97 + i)
        ax.set_title(f"({panel}) {tissue}", fontsize=9, fontweight="bold")
        # mark Bonferroni-significant cells
        ys, xs = np.where(grid.values < bonf)
        ax.scatter(xs, ys, marker="*", c="cyan", s=22, linewidths=0)

    cbar = fig.colorbar(im, ax=axes, shrink=0.6, aspect=30,
                        label=r"$-\log_{10}\,P$ (PC $\sim$ covariate)")
    cbar.ax.tick_params(labelsize=7)
    fig.suptitle(
        "STARNET expression PCs vs sample covariates "
        r"(cyan $\star$: Bonferroni-significant)",
        fontsize=10)
    for ext in ("png", "pdf"):
        fig.savefig(fig_dir / f"fig_covariate_check.{ext}")
    plt.close(fig)

    # ---- Figure 2 (fig_covariate_pc_scatter): PC1-PC2 coloured by Age, Sex and BMI ----
    # Written per tissue to the QC directory; the representative tissue is additionally
    # saved as the supplementary figure fig_covariate_pc_scatter.
    key_covs = [c for c in ("Age", "Sex", "BMI") if c in cov_names]
    pc1_var = {t: var_for_pc(long, t, "PC1") for t in tissues}
    pc2_var = {t: var_for_pc(long, t, "PC2") for t in tissues}
    for tissue in tissues:
        pcs = pcs_by_tissue[tissue]
        ph_t = aligned_ph[tissue]
        fig, axs = plt.subplots(1, len(key_covs),
                                figsize=(2.9 * len(key_covs), 2.9),
                                squeeze=False, constrained_layout=True)
        xlab = f"PC1 ({pc1_var[tissue] * 100:.1f}%)"
        ylab = f"PC2 ({pc2_var[tissue] * 100:.1f}%)"
        for j, cov in enumerate(key_covs):
            ax = axs[0][j]
            c = ph_t[cov].values.astype(float)
            if cov == "Sex":
                # binary covariate: two-colour discrete legend, not a bar
                for val, colour, label in (
                        (0.0, "#d55e00", "female"), (1.0, "#0072b2", "male")):
                    m = c == val
                    ax.scatter(pcs["PC1"][m], pcs["PC2"][m], c=colour, s=12,
                               edgecolors="none", label=label)
                ax.legend(title="Sex", fontsize=7, title_fontsize=7,
                          frameon=False, loc="best")
            else:
                sc = ax.scatter(pcs["PC1"], pcs["PC2"], c=c, cmap="viridis",
                                s=12, edgecolors="none")
                fig.colorbar(sc, ax=ax, shrink=0.85, label=COV_LABELS.get(cov, cov))
            ax.set_xlabel(xlab)
            ax.set_ylabel(ylab if j == 0 else "")
            ax.set_title(f"({chr(97 + j)}) {COV_LABELS.get(cov, cov)}", fontsize=9)
        fig.suptitle(f"{tissue}: leading expression PCs vs covariates",
                     fontsize=10)
        for ext in ("png", "pdf"):
            fig.savefig(out_dir / f"fig_pc_scatter_{tissue}.{ext}")
        if tissue == args.rep_tissue:
            for ext in ("png", "pdf"):
                fig.savefig(fig_dir / f"fig_covariate_pc_scatter.{ext}")
        plt.close(fig)

    # ---- SUMMARY.md ----
    sig = long[long["p_value"] < bonf].copy()
    sig = sig.sort_values("p_value")
    # focus on the covariates the `lm` step claims to remove
    age_sex = sig[sig["covariate"].isin(["Age", "Sex"])]
    lines = []
    lines.append("# STARNET expression covariate-adjustment check\n")
    lines.append(
        "**Question:** are the STARNET expression matrices "
        "(`*.EDAseq.gc.lm.or.RQN`) already covariate-adjusted, as the `lm` "
        "(linear-model covariate removal) step in the filename claims?\n")
    lines.append(
        "**Method:** for each tissue, compute the top "
        f"{args.n_pcs} principal components of the samples x genes expression "
        "matrix, then test each PC for association (Pearson r) with sample "
        "covariates from the STARNET cohort workbook. If the `lm` step removed "
        "age and sex, the leading PCs should show **no** significant "
        "association with Age/Sex.\n")
    lines.append(f"- tissues tested: {', '.join(tissues)}")
    lines.append(f"- covariates: {', '.join(cov_names)}")
    lines.append(f"- total tests: {n_tests}; Bonferroni alpha = {bonf:.2e}\n")
    lines.append("## Age / Sex associations (the covariates `lm` claims to remove)\n")
    if age_sex.empty:
        lines.append(
            "**None.** No PC in any tissue is associated with Age or Sex at "
            "Bonferroni significance -> consistent with the claim that age and "
            "sex were regressed out upstream.\n")
    else:
        lines.append(
            f"{len(age_sex)} Age/Sex associations survive Bonferroni "
            "(see table). Inspect whether these are leading PCs (high "
            "var_explained) — strong leading-PC association would contradict "
            "the adjustment claim.\n")
        lines.append(age_sex.head(20).to_markdown(index=False))
        lines.append("")
    lines.append("## All Bonferroni-significant PC-covariate associations\n")
    if sig.empty:
        lines.append("None.\n")
    else:
        lines.append(sig.head(40).to_markdown(index=False))
        lines.append("")
    (out_dir / "SUMMARY.md").write_text("\n".join(lines))
    print(f"[cov-test] wrote QC tables to {out_dir}/ and figures to {fig_dir}/")
    print(f"[cov-test] Age/Sex Bonferroni-significant hits: {len(age_sex)}")


if __name__ == "__main__":
    main()
