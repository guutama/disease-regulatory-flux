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

Inputs
  --weights    <tissue>.expr_model_weights.tsv.gz (gene_id, method, config, channel,
               variant_id, weight) -- one block per predictable gene
  --genotype   <tissue>.genotype_012.tsv.gz (variant_id, chromosome, position, ref, alt,
               then samples) -- supplies sd_j and sigma_g
  --gwas       <trait>.gwas.tsv.gz (variant_id, z, ...) -- ALT-aligned GWAS from Stage 8
  --tissue     --trait     labels written into the outputs
  --outdir

Output: <outdir>/association_<tissue>_<trait>.tsv, one row per gene
  gene_id  config  n_snp  n_cis  n_trans  n_gwas  z_twas  z_cis  z_trans
  sigma_g  p_value  p_adj  tissue  trait
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


def load_weights(path, delim, gene_col, variant_col):
    """gene_id -> (config, [(channel, variant_id, weight), ...]) in first-seen gene order."""
    genes: "OrderedDict[str, list]" = OrderedDict()
    config: dict[str, str] = {}
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
            genes.setdefault(g, []).append((f[chi], f[vi], float(f[wi])))
            config.setdefault(g, f[ci])
    return genes, config


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
    p.add_argument("--genotype", required=True, help="input: genotype_012 hard-call matrix")
    p.add_argument("--gwas", required=True, help="input: ALT-aligned GWAS (Stage 8)")
    p.add_argument("--tissue", required=True, help="tissue label")
    p.add_argument("--trait", required=True, help="trait label")
    p.add_argument("--outdir", required=True, help="output directory")
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

    genes, config = load_weights(args.weights, d, args.gene_col, args.variant_col)
    wanted = {v for rows in genes.values() for _, v, _ in rows}
    geno = load_genotype(args.genotype, d, args.variant_col, wanted)
    sd = {v: float(g.std()) for v, g in geno.items()}
    gwas_z = load_gwas_z(args.gwas, d, args.variant_col)

    records = []
    for gene, rows in genes.items():
        snps = [(ch, v, w) for ch, v, w in rows if v in geno]
        if not snps:
            continue
        pred = np.zeros_like(next(iter(geno.values())))
        for _, v, w in snps:
            pred = pred + w * geno[v]
        sigma_g = float(pred.std())
        if sigma_g == 0.0:
            continue
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
        records.append(dict(gene_id=gene, config=config[gene], n_snp=len(snps),
                            n_cis=n_cis, n_trans=n_trans, n_gwas=n_gwas,
                            z_twas=z_twas, z_cis=z_cis, z_trans=z_trans,
                            sigma_g=sigma_g, p_value=p_value))

    p_adj = bh_fdr([r["p_value"] for r in records])

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"association_{args.tissue}_{args.trait}.tsv")
    cols = ["gene_id", "config", "n_snp", "n_cis", "n_trans", "n_gwas",
            "z_twas", "z_cis", "z_trans", "sigma_g", "p_value", "p_adj", "tissue", "trait"]
    with open(out, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r, padj in zip(records, p_adj):
            r["p_adj"] = float(padj)
            r["tissue"], r["trait"] = args.tissue, args.trait
            fh.write("\t".join(
                f"{r[c]:.10g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")

    n_sig = int((p_adj < 0.05).sum())
    print(f"[twas_association] tissue={args.tissue} trait={args.trait} "
          f"tested={len(records)} sig(FDR<0.05)={n_sig}")


if __name__ == "__main__":
    main()
