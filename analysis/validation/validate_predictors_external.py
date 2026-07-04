#!/usr/bin/env python3
"""External validation of the expression predictors' LOO-R2, using the published weights only.

Our leave-one-out R2 (Supplementary Table S1) is an in-cohort claim about how well each gene's
channel-aware predictor reconstructs its expression. Anyone with genotypes and measured expression
from a matching tissue can test that claim in their own cohort WITHOUT retraining: the predictor is
a fixed linear function of dosage, so validating it is a matrix product plus a correlation.

For each gene, the predicted expression of sample i is

    yhat_i = sum_j  weight_j * dosage_ij ,

summed over ALL of the gene's Supplementary-Table-S9 SNPs -- both the cis channel and the trans
channel (the trans contribution is already baked into the per-SNP weights, so the external user
never needs the regulatory network). dosage_ij is the ALT-allele dosage (0-2) of SNP j, with the
ALT allele defined by S9's ref/alt.

The external readout is the squared Pearson correlation between predicted and measured expression,

    R2_external = cor(yhat, y_measured) ** 2 ,

which is invariant to the overall scale and offset of expression (a raw 1 - SSE/SST would be
distorted by cross-cohort unit differences). R2_external is compared to our reported loo_r2.

Inputs
    --weights      Supplementary Table S9 (per-SNP predictor weights with ref/alt), CSV[.gz]
    --loo          Supplementary Table S1 (supplies loo_r2 / predictable per gene-tissue), CSV
    --dosage       external dosage matrix: one row per variant, one column per sample, with an
                   id column (--dosage-id-col) and, to harmonise alleles, optional ref/alt columns
                   (--dosage-ref-col / --dosage-alt-col). Values are ALT-allele dosages (0-2).
    --expression   external expression matrix: one row per gene (Ensembl id), one column per
                   sample. Should be preprocessed like the discovery expression (covariate-adjusted
                   and normalised) so it lives in the space the predictor targets.
    --tissue       tissue whose S9 weights to apply (must match the external tissue).

Notes on alignment (the external user is responsible for these):
    - Build must match S9 (GRCh37). Variants are matched by the id column.
    - If ref/alt columns are given, dosage is flipped (2 - dosage) when the external ALT is S9's
      REF, and variants with non-matching alleles are dropped. Without them the dosage is assumed
      already ALT-oriented to S9.
    - Genes are matched by version-stripped Ensembl id; SNPs absent externally are dropped and
      per-sample missing dosages are mean-imputed.

Output: a per-gene table (gene, gene_id, tissue, config, n_snp_total, n_snp_used, n_snp_missing,
n_samples, pearson_r, r2_external, loo_r2, predictable) and, with --fig, a loo_r2-vs-R2_external
scatter in the shared figures directory.
"""
import argparse
import gzip
import os
from pathlib import Path

import numpy as np
import pandas as pd

# The scatter figure goes to the shared figures directory: results/figures/ by default.
FIGDIR = Path(os.environ.get("FLUX_FIGURES", Path(__file__).resolve().parents[2] / "results" / "figures"))


def base_id(g):
    """Version-stripped Ensembl gene id."""
    return str(g).split(".")[0]


def load_loo(path, tissue):
    """gene_id -> (loo_r2, predictable) from Supplementary Table S1, for one tissue."""
    df = pd.read_csv(path, usecols=["gene_id", "tissue", "best_loo_r2", "predictable"])
    df = df[df["tissue"] == tissue]
    pred = df["predictable"].astype(str).str.lower() == "true"
    return {gid: (r2, bool(p)) for gid, r2, p in zip(df["gene_id"], df["best_loo_r2"], pred)}


def read_dosage(path, wanted_ids, id_col, ref_col=None, alt_col=None):
    """Stream the external dosage matrix, keeping only the wanted variant ids. Returns
    (samples, raw) where raw maps id -> (ext_ref, ext_alt, dosage array over samples). ext_ref /
    ext_alt are None when no allele columns are given."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        meta = [id_col] + [c for c in (ref_col, alt_col) if c]
        for c in meta:
            if c not in idx:
                raise SystemExit(f"dosage matrix has no column '{c}'")
        i_id = idx[id_col]
        i_ref = idx[ref_col] if ref_col else None
        i_alt = idx[alt_col] if alt_col else None
        sample_cols = [i for i, c in enumerate(header) if c not in meta]
        samples = [header[i] for i in sample_cols]
        raw = {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            vid = f[i_id]
            if vid not in wanted_ids:
                continue
            dos = np.array([np.nan if f[i] in ("", "NA", "nan") else float(f[i])
                            for i in sample_cols], dtype=float)
            er = f[i_ref] if i_ref is not None else None
            ea = f[i_alt] if i_alt is not None else None
            raw[vid] = (er, ea, dos)
    return samples, raw


def harmonise(raw, alleles):
    """Orient external dosages to the S9 ALT allele. alleles maps id -> (s9_ref, s9_alt). Returns
    (oriented, n_flipped, n_mismatch): oriented maps id -> dosage array on the S9 ALT allele."""
    oriented, n_flipped, n_mismatch = {}, 0, 0
    for vid, (er, ea, dos) in raw.items():
        sr, sa = alleles[vid]
        if er is None or ea is None:                       # no allele info: assume aligned
            oriented[vid] = dos
        elif (er, ea) == (sr, sa):
            oriented[vid] = dos
        elif (er, ea) == (sa, sr):
            oriented[vid] = 2.0 - dos                       # external counts the S9 REF
            n_flipped += 1
        else:
            n_mismatch += 1                                 # different variant / strand: drop
    return oriented, n_flipped, n_mismatch


def _impute(mat):
    """Replace per-SNP missing dosages with that SNP's mean; drop SNPs that are entirely missing."""
    keep = ~np.all(np.isnan(mat), axis=1)
    mat = mat[keep]
    for r in range(mat.shape[0]):
        row = mat[r]
        m = np.isnan(row)
        if m.any():
            row[m] = np.nanmean(row)
    return mat, keep


def compute(weights, dosage, expr, loo, min_snp=1):
    """Score every gene in `weights` (S9 rows for one tissue). `dosage` is a DataFrame indexed by
    variant id with one column per sample (ALT-oriented); `expr` is a DataFrame indexed by
    version-stripped Ensembl id with the same sample columns; `loo` maps gene_id -> (loo_r2,
    predictable). Returns one row per gene with the external squared correlation."""
    dos_ids = set(dosage.index)
    rows = []
    for gene_id, g in weights.groupby("gene_id", sort=False):
        base = base_id(gene_id)
        loo_r2, predictable = loo.get(gene_id, (np.nan, False))
        rec = {"gene": g["gene"].iloc[0], "gene_id": gene_id, "tissue": g["tissue"].iloc[0],
               "config": g["config"].iloc[0], "n_snp_total": len(g),
               "n_snp_used": 0, "n_snp_missing": len(g), "n_samples": 0,
               "pearson_r": np.nan, "r2_external": np.nan,
               "loo_r2": loo_r2, "predictable": predictable}
        present = g[g["rs_id"].isin(dos_ids)]
        rec["n_snp_used"] = len(present)
        rec["n_snp_missing"] = len(g) - len(present)
        if base not in expr.index or len(present) < min_snp:
            rows.append(rec)
            continue
        mat = dosage.loc[present["rs_id"]].to_numpy(dtype=float)      # (n_snp, n_sample)
        w = present["weight"].to_numpy(dtype=float)
        mat, keep = _impute(mat)
        w = w[keep]
        if mat.shape[0] < min_snp:
            rows.append(rec)
            continue
        yhat = w @ mat                                                # (n_sample,)
        y = expr.loc[base].to_numpy(dtype=float)
        mask = ~np.isnan(y) & ~np.isnan(yhat)
        rec["n_samples"] = int(mask.sum())
        if mask.sum() >= 3 and np.nanstd(yhat[mask]) > 0 and np.nanstd(y[mask]) > 0:
            r = float(np.corrcoef(yhat[mask], y[mask])[0, 1])
            rec["pearson_r"] = r
            rec["r2_external"] = r * r
        rows.append(rec)
    cols = ["gene", "gene_id", "tissue", "config", "n_snp_total", "n_snp_used", "n_snp_missing",
            "n_samples", "pearson_r", "r2_external", "loo_r2", "predictable"]
    return pd.DataFrame(rows, columns=cols)


def scatter(out, fig_path):
    """loo_r2 (x) vs R2_external (y) scatter over scored genes, highlighting predictable genes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = out.dropna(subset=["r2_external", "loo_r2"])
    if d.empty:
        return
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 9, "savefig.dpi": 300})
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    other = d[~d["predictable"].astype(bool)]
    pred = d[d["predictable"].astype(bool)]
    ax.scatter(other["loo_r2"], other["r2_external"], s=10, c="#c2c5cc",
               edgecolors="none", label="not predictable")
    ax.scatter(pred["loo_r2"], pred["r2_external"], s=12, c="#2C6FB5",
               edgecolors="none", label="predictable (loo $R^2\\geq$0.01)")
    lim = float(np.nanmax([d["loo_r2"].max(), d["r2_external"].max()])) * 1.05
    ax.plot([0, lim], [0, lim], ls="--", lw=0.8, color="0.5")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
    ax.set_xlabel("discovery LOO $R^2$ (S1)")
    ax.set_ylabel("external $R^2$  (cor$^2$)")
    if len(pred) >= 2:
        rr = np.corrcoef(pred["loo_r2"], pred["r2_external"])[0, 1]
        ax.text(0.04, 0.96, f"$r={rr:.2f}$\n$n={len(pred)}$", transform=ax.transAxes,
                va="top", ha="left", fontsize=8)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(str(fig_path).replace(".pdf", f".{ext}"), bbox_inches="tight")
    plt.close(fig)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--weights", required=True, help="Supplementary Table S9 (CSV[.gz])")
    ap.add_argument("--loo", required=True, help="Supplementary Table S1 (CSV)")
    ap.add_argument("--dosage", required=True, help="external dosage matrix (TSV[.gz])")
    ap.add_argument("--expression", required=True, help="external expression matrix (TSV[.gz])")
    ap.add_argument("--tissue", required=True, help="tissue whose S9 weights to apply")
    ap.add_argument("--out", required=True, help="output per-gene validation table (CSV)")
    ap.add_argument("--dosage-id-col", default="rs_id")
    ap.add_argument("--dosage-ref-col", help="external REF-allele column (enables harmonisation)")
    ap.add_argument("--dosage-alt-col", help="external ALT-allele column (enables harmonisation)")
    ap.add_argument("--expr-id-col", default="gene_id")
    ap.add_argument("--min-snp", type=int, default=1,
                    help="minimum external SNPs to score a gene (default: 1)")
    ap.add_argument("--fig", nargs="?", const=str(FIGDIR / "fig_external_validation.pdf"),
                    help="write the loo-vs-external scatter (default path in the shared "
                         "figures directory)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    weights = pd.read_csv(args.weights, usecols=["gene", "gene_id", "tissue", "config",
                                                 "rs_id", "ref", "alt", "weight"])
    weights = weights[weights["tissue"] == args.tissue].copy()
    if weights.empty:
        raise SystemExit(f"no S9 weights for tissue '{args.tissue}'")
    alleles = {r.rs_id: (r.ref, r.alt) for r in
               weights[["rs_id", "ref", "alt"]].drop_duplicates().itertuples(index=False)}
    wanted = set(weights["rs_id"])

    samples_d, raw = read_dosage(args.dosage, wanted, args.dosage_id_col,
                                 args.dosage_ref_col, args.dosage_alt_col)
    oriented, n_flip, n_mm = harmonise(raw, alleles)
    print(f"[ext] {len(oriented)}/{len(wanted)} predictor SNPs found in the external dosage "
          f"({n_flip} allele-flipped, {n_mm} allele-mismatched and dropped)")

    expr = pd.read_csv(args.expression, sep=None, engine="python")
    eid = args.expr_id_col
    expr[eid] = expr[eid].map(base_id)
    expr = expr.drop_duplicates(subset=[eid]).set_index(eid)
    samples = [s for s in samples_d if s in expr.columns]
    if len(samples) < 3:
        raise SystemExit(f"only {len(samples)} shared samples between dosage and expression")
    print(f"[ext] {len(samples)} samples shared between dosage and expression")

    dos_df = pd.DataFrame({vid: pd.Series(dict(zip(samples_d, arr))) for vid, arr in oriented.items()}).T
    dos_df = dos_df[samples]
    expr = expr[samples].apply(pd.to_numeric, errors="coerce")

    loo = load_loo(args.loo, args.tissue)
    out = compute(weights, dos_df, expr, loo, min_snp=args.min_snp)
    out.to_csv(args.out, index=False)

    scored = out.dropna(subset=["r2_external"])
    pred = scored[scored["predictable"].astype(bool)]
    print(f"[ext] scored {len(scored)} genes; mean external R2 = "
          f"{scored['r2_external'].mean():.4f} (all), {pred['r2_external'].mean():.4f} (predictable)")
    if len(pred) >= 2:
        keep = pred.dropna(subset=["loo_r2"])
        r = np.corrcoef(keep["loo_r2"], keep["r2_external"])[0, 1] if len(keep) >= 2 else np.nan
        frac = float((pred["r2_external"] >= 0.01).mean())
        print(f"[ext] predictable genes: cor(loo_r2, external_R2) = {r:.3f}; "
              f"{100*frac:.0f}% keep external R2 >= 0.01")
    print(f"[ext] wrote {args.out}")

    if args.fig:
        scatter(out, Path(args.fig))
        print(f"[ext] wrote {args.fig}")


if __name__ == "__main__":
    main()
