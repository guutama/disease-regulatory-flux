#!/usr/bin/env python3
"""Stage 9 -- summary-statistic TWAS association with a cis/trans decomposition.

Each gene's selected expression predictor (Stage 7) carries ALT-dosage weights on its cis
and trans SNPs. Combined with an ALT-aligned GWAS (Stage 8), the per-gene association is

    z_TWAS = ( sum_j  w_j * sd_j * z_gwas_j ) / sigma_g ,

where w_j is the raw ALT-dosage weight of SNP j, sd_j its genotype standard deviation,
z_gwas_j the GWAS z of its ALT allele, and sigma_g = sd( G w ) the spread of the predicted
expression over the genotypes (the LD/variance term -- no external panel needed). Restricting
the numerator to one channel gives the cis and trans components, which sum to z_TWAS. The
two-sided p-value is from the standard normal; a Benjamini-Hochberg FDR is applied per tissue.

Predictor-stability diagnostics accompany every gene, because a small or degenerate predictor
can inflate |z_TWAS| off null SNPs. Two dimensionless quantities are reported:

    max_amp = max_j |w_j sd_j| / sigma_g       single-SNP leverage ("explosion")
    cancel  = sqrt( sum_j (w_j sd_j)^2 ) / sigma_g   weight cancellation ("instability")

max_amp is how much a single SNP's dosage-scaled weight exceeds the whole predictor's SD, so
that SNP alone drives z ~ max_amp * z_gwas_j when sigma_g is too small next to it. cancel is how
much smaller the realised sigma_g is than the independent-SNP scale, i.e. how much the SNPs
cancel in G w and shrink the denominator. A gene is flagged unstable when max_amp >= --amp-thr
or cancel >= --cancel-thr (both default 3; a single-SNP predictor has max_amp = cancel = 1 by
construction). Separately, genes whose predictor SD falls below --sigma-g-floor (an absolute
guard against a degenerate denominator) are dropped before FDR, alongside the sigma_g == 0 case.

What cancel really is, and what it does NOT claim. Writing r_j = w_j sd_j for the scaled
weights and R for the SNP LD correlation matrix, sigma_g^2 = r' R r and sum_j (w_j sd_j)^2 = r'r,
so cancel^2 = (r'r) / (r' R r) = 1 / Rayleigh(R, r): cancel is the reciprocal Rayleigh quotient
of the LD matrix along the weight direction -- a per-gene condition number of the burden
denominator w'Sigma w (the same quadratic form FUSION-style TWAS puts under the square root).
Because tr(R) = n, cancel < 1 means r lies along high-eigenvalue (correlated, reinforcing)
directions and cancel > 1 along low-eigenvalue (near-collinear, cancelling) ones, bounded by
cancel^2 <= 1/lambda_min(R). Crucially, cancel is a conditioning/sensitivity flag, not evidence
of a false positive: under the exactly-specified null z_gwas ~ N(0, R), z_TWAS = r'z / sigma_g
has variance 1 for ANY cancel -- sigma_g is the correct normalization. cancel is precisely the
factor by which any mismatch between the LD in sigma_g and the true joint law of the GWAS z's is
amplified. Here sigma_g uses the in-sample genotypes the predictor was fit on, so there is no
reference-panel mismatch on the denominator; the only exposure of a high-cancel gene is between
that in-sample LD and the external GWAS cohort's LD (plus pleiotropy / real multi-SNP effects).
That is why max_amp/cancel/unstable are reported as columns rather than row-dropped: the hard
filter that removes flagged genes from downstream results is a conservative robustness policy
(we cannot certify per-gene in-sample-vs-GWAS LD agreement), not a verdict of spuriousness.

Inputs
  --weights    <tissue>.expr_model_weights.tsv.gz (gene_id, method, config, channel,
               variant_id, weight) -- one block per predictable gene
  --genotype   <tissue>.genotype_012.tsv.gz (variant_id, chromosome, position, ref, alt,
               then samples) -- supplies sd_j and sigma_g
  --gwas       <trait>.gwas.tsv.gz (variant_id, z, ...) -- ALT-aligned GWAS from Stage 8
  --tissue     --trait     labels written into the outputs
  --outdir

A trans-channel significance floor (--trans-sigma-floor) reflects that a trans predictor with
a small sigma_g inflates |z_trans| off the small denominator: a gene whose winning predictor
uses the trans channel (n_trans > 0) and whose sigma_g is at or below the floor is not called
significant (its p_adj is forced to 1 and trans_ok is False), because the trans association is
not trusted below the floor. The raw two-sided p_value is left untouched.

Output: <outdir>/association_<tissue>_<trait>.tsv, one row per gene
  gene_id  config  n_snp  n_cis  n_trans  n_gwas  z_twas  z_cis  z_trans
  sigma_g  max_amp  cancel  unstable  p_value  p_adj  trans_ok  tissue  trait
"""
from __future__ import annotations

import argparse
import gzip
import math
import os
from collections import OrderedDict

import numpy as np

_SQRT2 = math.sqrt(2.0)


def _open(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _split(line: str, delim: str) -> list[str]:
    return line.rstrip("\n").split(delim)


def _col(header: list[str], name: str, what: str) -> int:
    if name not in header:
        raise SystemExit(f"The {what} file must have a '{name}' column; it has {header}.")
    return header.index(name)


def load_selected(path, delim, gene_col):
    """gene_id -> selected (winning) config, for the predictable genes only."""
    sel: dict[str, str] = {}
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        gi = _col(header, gene_col, "selected")
        ci = _col(header, "best_config", "selected")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            sel[f[gi]] = f[ci]
    return sel


def load_weights(path, delim, gene_col, variant_col, selected):
    """gene_id -> [(channel, variant_id, weight), ...] in first-seen gene order.

    The weights file holds every fitted configuration; we keep only the rows of each
    predictable gene's selected config, so the predictor is a single coherent model."""
    genes: "OrderedDict[str, list]" = OrderedDict()
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        gi = _col(header, gene_col, "weights")
        ci = _col(header, "config", "weights")
        chi = _col(header, "channel", "weights")
        vi = _col(header, variant_col, "weights")
        wi = _col(header, "weight", "weights")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            g = f[gi]
            if selected.get(g) != f[ci]:
                continue
            genes.setdefault(g, []).append((f[chi], f[vi], float(f[wi])))
    return genes


def load_genotype(path, delim, variant_col, wanted):
    """variant_id -> dosage vector (float), for the wanted variants. Genotype columns:
    variant_id, chromosome, position, ref, alt, then one per sample."""
    geno: dict[str, np.ndarray] = {}
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        if header[:5] != [variant_col, "chromosome", "position", "ref", "alt"]:
            raise SystemExit(
                f"The genotype matrix must start with the columns '{variant_col}', "
                f"'chromosome', 'position', 'ref', 'alt', but it starts with {header[:5]}.")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            if f[0] in wanted and f[0] not in geno:
                geno[f[0]] = np.array(f[5:], dtype=np.float64)
    return geno


def load_gwas_z(path, delim, variant_col):
    """variant_id -> ALT-aligned z (float)."""
    z: dict[str, float] = {}
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        vi = _col(header, variant_col, "GWAS")
        zi = _col(header, "z", "GWAS")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            try:
                z[f[vi]] = float(f[zi])
            except ValueError:
                continue
        return z


def bh_fdr(pvals: list[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(pvals, dtype=np.float64)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0.0, 1.0)
    return adj


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--weights", required=True, help="input: expr_model_weights matrix")
    p.add_argument("--selected", required=True,
                   help="input: expr_model_selected table (each predictable gene's best config)")
    p.add_argument("--genotype", required=True, help="input: genotype_012 hard-call matrix")
    p.add_argument("--gwas", required=True, help="input: ALT-aligned GWAS (Stage 8)")
    p.add_argument("--tissue", required=True, help="tissue label")
    p.add_argument("--trait", required=True, help="trait label")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--amp-thr", type=float, default=3.0,
                   help="unstable if max_amp = max_j |w_j sd_j|/sigma_g >= this (default: 3)")
    p.add_argument("--cancel-thr", type=float, default=3.0,
                   help="unstable if cancel = sqrt(sum_j (w_j sd_j)^2)/sigma_g >= this (default: 3)")
    p.add_argument("--sigma-g-floor", type=float, default=0.0,
                   help="drop genes whose predictor SD sigma_g is below this absolute floor "
                        "(degenerate-denominator guard; default: 0 = only drop sigma_g == 0)")
    p.add_argument("--trans-sigma-floor", type=float, default=0.0,
                   help="a trans-using gene (n_trans > 0) with sigma_g <= this is not called "
                        "significant (p_adj forced to 1, trans_ok=False); distrusts the trans "
                        "channel below the floor (default: 0 = off)")
    p.add_argument("--delimiter", default="\t", help="field delimiter (default: tab)")
    p.add_argument("--gene-col", default="gene_id", help="gene-id column (default: gene_id)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column (default: variant_id)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    selected = load_selected(args.selected, d, args.gene_col)
    genes = load_weights(args.weights, d, args.gene_col, args.variant_col, selected)
    wanted = {v for rows in genes.values() for _, v, _ in rows}
    geno = load_genotype(args.genotype, d, args.variant_col, wanted)
    sd = {v: float(g.std()) for v, g in geno.items()}
    gwas_z = load_gwas_z(args.gwas, d, args.variant_col)

    records = []
    n_degenerate = 0
    for gene, rows in genes.items():
        snps = [(ch, v, w) for ch, v, w in rows if v in geno]
        if not snps:
            continue
        pred = np.zeros_like(next(iter(geno.values())))
        for _, v, w in snps:
            pred = pred + w * geno[v]
        sigma_g = float(pred.std())
        if sigma_g == 0.0 or sigma_g < args.sigma_g_floor:
            n_degenerate += 1          # degenerate denominator -> drop before FDR
            continue
        # predictor-stability diagnostics (GWAS-independent properties of the weights)
        wsd = np.array([w * sd[v] for _, v, w in snps])
        max_amp = float(np.max(np.abs(wsd)) / sigma_g)
        cancel = float(np.sqrt(np.sum(wsd ** 2)) / sigma_g)
        unstable = (max_amp >= args.amp_thr) or (cancel >= args.cancel_thr)
        num_cis = num_trans = 0.0
        n_cis = n_trans = n_gwas = 0
        for ch, v, w in snps:
            if ch == "cis":
                n_cis += 1
            elif ch == "trans":
                n_trans += 1
            if v in gwas_z:
                contrib = w * sd[v] * gwas_z[v]
                if ch == "cis":
                    num_cis += contrib
                elif ch == "trans":
                    num_trans += contrib
                n_gwas += 1
        if n_gwas == 0:
            continue
        z_cis = num_cis / sigma_g
        z_trans = num_trans / sigma_g
        z_twas = z_cis + z_trans
        p_value = math.erfc(abs(z_twas) / _SQRT2)
        records.append(dict(gene_id=gene, config=selected[gene], n_snp=len(snps),
                            n_cis=n_cis, n_trans=n_trans, n_gwas=n_gwas,
                            z_twas=z_twas, z_cis=z_cis, z_trans=z_trans,
                            sigma_g=sigma_g, max_amp=max_amp, cancel=cancel,
                            unstable=unstable, p_value=p_value))

    p_adj = bh_fdr([r["p_value"] for r in records])
    # trans-channel significance floor: a distrusted small-sigma_g trans predictor is excluded
    # from the significant set (p_adj -> 1); trans_ok records whether the gene passed the rule.
    for i, r in enumerate(records):
        r["trans_ok"] = not (r["n_trans"] > 0 and r["sigma_g"] <= args.trans_sigma_floor)
        if not r["trans_ok"]:
            p_adj[i] = 1.0

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"association_{args.tissue}_{args.trait}.tsv")
    cols = ["gene_id", "config", "n_snp", "n_cis", "n_trans", "n_gwas",
            "z_twas", "z_cis", "z_trans", "sigma_g", "max_amp", "cancel", "unstable",
            "p_value", "p_adj", "trans_ok", "tissue", "trait"]
    with open(out, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r, padj in zip(records, p_adj):
            r["p_adj"] = float(padj)
            r["tissue"], r["trait"] = args.tissue, args.trait
            fh.write("\t".join(
                f"{r[c]:.10g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

    n_sig = int((p_adj < 0.05).sum())
    n_unstable = sum(r["unstable"] for r in records)
    n_trans_floored = sum(1 for r in records if not r["trans_ok"])
    print(f"[twas_association] tissue={args.tissue} trait={args.trait} "
          f"tested={len(records)} sig(FDR<0.05)={n_sig} unstable={n_unstable} "
          f"dropped_low_sigma={n_degenerate} trans_floored={n_trans_floored}")


if __name__ == "__main__":
    main()
