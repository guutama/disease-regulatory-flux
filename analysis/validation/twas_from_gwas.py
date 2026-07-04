#!/usr/bin/env python3
"""External summary-statistic TWAS from the published predictor weights.

Applies the disease-regulatory-flux expression predictors (Supplementary Table S9) to any
GWAS to obtain a per-gene, per-tissue TWAS statistic, decomposed into its cis and trans
parts. Only published, non-restricted tables are needed, so anyone can validate or extend
the associations on a new trait without access to the individual-level STARNET data.

For each gene's selected predictor, the association is

    z_TWAS = ( sum_j  weight_j * sd_j * z_gwas_j ) / sigma_g ,

where weight_j and sd_j are the SNP's weight and training-sample ALT-dosage SD (both from
S9), z_gwas_j is the GWAS z of its ALT allele, and sigma_g is the SD of the predicted
expression (from Supplementary Table S2). Restricting the sum to one channel gives the cis
and trans components, which add up to z_TWAS. SNPs whose sd is blank in S9 (outside the
association's genotype panel) or absent from the GWAS are skipped, exactly as in discovery.
Two-sided p-values come from the standard normal; a Benjamini-Hochberg FDR is applied per
tissue.

The GWAS must first be aligned to the S9 ALT alleles. Do this once with the pipeline's
harmoniser, using S9 as the variant reference:

    bin/harmonise_gwas.py --gwas <raw_gwas> \\
        --variants Supplementary_Table_S9_predictor_weights.csv.gz \\
        --variant-col rs_id --out <trait>.aligned.tsv.gz

That leaves an rs_id + z table already on the ALT allele; no access to the restricted
STARNET individual-level data is required.

Usage:
  python analysis/validation/twas_from_gwas.py \\
      --weights   Supplementary_Table_S9_predictor_weights.csv.gz \\
      --sigma-g   Supplementary_Table_S2_twas_association.csv \\
      --gwas      <trait>.aligned.tsv.gz \\
      --out       <trait>.external_twas.csv \\
      [--discovery Supplementary_Table_S2_twas_association.csv]

With --discovery, the discovery z_TWAS is joined on and a sign-concordance flag is added,
so the external trait can be compared against the reported CAD associations directly.
"""
import argparse
import math

import pandas as pd

_SQRT2 = math.sqrt(2.0)

OUT_COLS = [
    "gene", "gene_id", "tissue", "config",
    "n_snp", "n_cis", "n_trans", "n_gwas",
    "z_twas", "z_cis", "z_trans", "sigma_g", "p", "p_adj",
]


def two_sided_p(z):
    """Two-sided standard-normal p-value for a z-score."""
    return math.erfc(abs(z) / _SQRT2)


def bh_fdr(p):
    """Benjamini-Hochberg adjusted p-values for a Series, returned aligned to its index."""
    s = p.sort_values()
    n = len(s)
    ranked = (s.to_numpy() * n / (pd.Series(range(1, n + 1))).to_numpy())
    # enforce monotonicity from the largest p downwards
    adj = pd.Series(ranked, index=s.index)[::-1].cummin()[::-1].clip(upper=1.0)
    return adj.reindex(p.index)


def load_sigma_g(path):
    """(gene_id, tissue) -> sigma_g from the association table (Supplementary Table S2)."""
    df = pd.read_csv(path, usecols=["gene_id", "tissue", "sigma_g"])
    return df.set_index(["gene_id", "tissue"])["sigma_g"].to_dict()


def load_discovery(path):
    """(gene_id, tissue) -> discovery z_twas from Supplementary Table S2, for comparison."""
    df = pd.read_csv(path, usecols=["gene_id", "tissue", "z_twas"])
    return df.set_index(["gene_id", "tissue"])["z_twas"].to_dict()


def compute(weights, gwas_z, sigma_g):
    """Combine the per-SNP weights with the aligned GWAS into one row per gene-tissue."""
    w = weights.dropna(subset=["sd"]).copy()          # SNPs usable for summary TWAS
    w["z_gwas"] = w["rs_id"].map(gwas_z)
    w["contrib"] = w["weight"] * w["sd"] * w["z_gwas"]

    rows = []
    for (gene_id, tissue), g in w.groupby(["gene_id", "tissue"], sort=False):
        sg = sigma_g.get((gene_id, tissue))
        if sg is None or sg == 0 or pd.isna(sg):
            continue
        used = g.dropna(subset=["z_gwas"])
        if used.empty:
            continue
        z_cis = used.loc[used.channel == "cis", "contrib"].sum() / sg
        z_trans = used.loc[used.channel == "trans", "contrib"].sum() / sg
        z_twas = z_cis + z_trans
        rows.append({
            "gene": g["gene"].iloc[0],
            "gene_id": gene_id,
            "tissue": tissue,
            "config": g["config"].iloc[0],
            "n_snp": len(g),
            "n_cis": int((used.channel == "cis").sum()),
            "n_trans": int((used.channel == "trans").sum()),
            "n_gwas": len(used),
            "z_twas": z_twas,
            "z_cis": z_cis,
            "z_trans": z_trans,
            "sigma_g": sg,
            "p": two_sided_p(z_twas),
        })
    out = pd.DataFrame(rows, columns=[c for c in OUT_COLS if c != "p_adj"])
    if out.empty:
        return out.assign(p_adj=[])
    out["p_adj"] = out.groupby("tissue")["p"].transform(bh_fdr)
    return out[OUT_COLS]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--weights", required=True,
                    help="Supplementary Table S9 (per-SNP predictor weights with sd)")
    ap.add_argument("--sigma-g", required=True,
                    help="Supplementary Table S2 (supplies sigma_g per gene-tissue)")
    ap.add_argument("--gwas", required=True,
                    help="GWAS aligned to the S9 ALT alleles (rs_id + z columns)")
    ap.add_argument("--out", required=True, help="output per-gene TWAS table (CSV)")
    ap.add_argument("--discovery",
                    help="optional Supplementary Table S2 to join the discovery z_twas and "
                         "flag sign concordance")
    ap.add_argument("--gwas-variant-col", default="rs_id",
                    help="variant-id column in the GWAS (default: rs_id)")
    ap.add_argument("--gwas-z-col", default="z",
                    help="z-score column in the GWAS (default: z)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    weights = pd.read_csv(args.weights)
    gwas = pd.read_csv(args.gwas, sep=None, engine="python",
                       usecols=[args.gwas_variant_col, args.gwas_z_col])
    gwas_z = dict(zip(gwas[args.gwas_variant_col], gwas[args.gwas_z_col]))
    sigma_g = load_sigma_g(args.sigma_g)

    out = compute(weights, gwas_z, sigma_g)

    if args.discovery:
        disc = load_discovery(args.discovery)
        out["z_discovery"] = [disc.get((gid, t)) for gid, t in zip(out.gene_id, out.tissue)]
        out["concordant"] = ((out["z_twas"] > 0) == (out["z_discovery"] > 0)) \
            & out["z_discovery"].notna()

    out.to_csv(args.out, index=False)
    n_sig = int((out["p_adj"] < 0.05).sum()) if len(out) else 0
    print(f"external TWAS: {len(out)} gene-tissue rows, {n_sig} significant (FDR<0.05) "
          f"-> {args.out}")


if __name__ == "__main__":
    main()
