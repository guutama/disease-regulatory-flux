#!/usr/bin/env python3
"""
LD-prune each gene's cis-eQTL SNP set for one tissue and write a hard-call genotype matrix
for the SNPs that survive.

Neighbouring SNPs are often in strong linkage disequilibrium (LD), so a gene's cis-eQTL set
carries redundant, highly correlated variants. This step thins each gene's set: it walks the
gene's cis SNPs from the most to the least significant eQTL and keeps a SNP only if it is not
in LD --- squared correlation above a threshold --- with a SNP already kept for that gene.
The result is a smaller, less redundant set of cis features per gene.

LD is measured as the squared Pearson correlation between SNPs on the tissue's own samples,
using genotypes rounded to hard calls (0, 1 or 2 copies of the alternate allele). A SNP that
is constant across the samples carries no information and is dropped.

Inputs:
  --cis-features  the per-gene cis-SNP set (the cis-eQTL feature-selection step): gene_id,
                  variant_id.
  --genotype      the harmonised dosage matrix (variant_id followed by one column per
                  sample); dosages are rounded to 0/1/2 here.
  --eqtl          the harmonised cis-eQTL table, used only to order each gene's SNPs by
                  ascending p-value before clumping; optional --- without a p-value column
                  the SNPs are clumped in the order they appear in the cis-feature set.

Output (in --outdir), for tissue <tissue>:
  <tissue>.cis_snps_pruned.tsv.gz   gene_id, variant_id (the SNPs kept per gene);
  <tissue>.genotype_012.tsv.gz      variant_id followed by the sample columns, 0/1/2 hard
                                    calls, restricted to the kept SNPs across all genes;
  <tissue>.ld_prune_summary.tsv     per gene: cis SNPs in, kept, and dropped for LD;
  <tissue>.ld_prune_qc.tsv          tissue-level totals.

Column names and the field delimiters are options, so a new dataset is supported by adjusting
configuration rather than editing code.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict
from typing import TextIO

import numpy as np


# ----------------------------------------------------------------------------- io helpers
def _open(path: str, mode: str = "rt") -> TextIO:
    """Open transparently, gzip-decompressing iff the name ends in .gz; '-' is std stream."""
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _split(line: str, delim: str) -> list[str]:
    return line.rstrip("\n").split(delim)


def _col(header: list[str], name: str, where: str) -> int:
    try:
        return header.index(name)
    except ValueError:
        raise SystemExit(f"The {where} file is missing the required column '{name}'.")


# ------------------------------------------------------------------------ LD greedy clump
def greedy_clump(order: list[str], rows: dict[str, np.ndarray], r2_threshold: float):
    """Greedily keep SNPs that are not in LD with an already-kept one.

    `order` is the variant ids to consider, most-significant first; `rows` maps a variant id
    to its standardised genotype vector (mean 0, unit norm) or None if it is constant. Walk
    `order`; keep a SNP unless its squared correlation with any kept SNP exceeds
    `r2_threshold`. Returns (kept, n_dropped_constant, n_dropped_ld)."""
    kept: list[str] = []
    kept_vecs: list[np.ndarray] = []
    n_const = n_ld = 0
    thr = r2_threshold ** 0.5            # compare |r| to sqrt(threshold) to avoid squaring
    for v in order:
        vec = rows.get(v)
        if vec is None:                  # constant SNP: no information, undefined correlation
            n_const += 1
            continue
        if any(abs(float(vec @ kv)) > thr for kv in kept_vecs):   # vec.kv = Pearson r (unit norm)
            n_ld += 1
            continue
        kept.append(v)
        kept_vecs.append(vec)
    return kept, n_const, n_ld


def standardise(dosage: np.ndarray):
    """Return the mean-centred, unit-norm genotype vector, or None if the SNP is constant."""
    x = dosage.astype(np.float64)
    x -= x.mean()
    norm = np.sqrt(x @ x)
    if norm == 0.0:
        return None
    return x / norm


# ----------------------------------------------------------------------------------- main
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cis-features", required=True,
                   help="input: per-gene cis-SNP set (gene_id, variant_id)")
    p.add_argument("--genotype", required=True,
                   help="input: harmonised dosage matrix "
                        "(variant_id, chromosome, position, ref, alt, then sample columns)")
    p.add_argument("--eqtl", default=None,
                   help="input: harmonised eQTL table, used to order SNPs by p-value "
                        "(optional)")
    p.add_argument("--tissue", required=True, help="tissue label, used to name the outputs")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--ld-r2-threshold", type=float, default=0.7,
                   help="keep a SNP unless its r^2 with an already-kept SNP exceeds this "
                        "(default: 0.7)")
    # format options
    p.add_argument("--delimiter", default="\t",
                   help="field delimiter of the input tables (default: tab)")
    p.add_argument("--gene-col", default="gene_id",
                   help="gene-id column in the cis-feature and eQTL tables (default: gene_id)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column in all input tables (default: variant_id)")
    p.add_argument("--p-col", default="pvalue",
                   help="p-value column in the eQTL table, if present (default: pvalue)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    # ---- per-gene cis SNP set, and the union universe ----------------------------------
    gene_snps: dict[str, list[str]] = defaultdict(list)
    universe: set[str] = set()
    with _open(args.cis_features) as fh:
        header = _split(fh.readline(), d)
        gi = _col(header, args.gene_col, "cis-feature")
        vi = _col(header, args.variant_col, "cis-feature")
        for line in fh:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            gene_snps[f[gi]].append(f[vi])
            universe.add(f[vi])
    if not universe:
        raise SystemExit("The cis-feature set is empty, so there is nothing to prune.")

    # ---- p-values for ordering (optional) ----------------------------------------------
    pval: dict[tuple[str, str], float] = {}
    if args.eqtl:
        with _open(args.eqtl) as fh:
            header = _split(fh.readline(), d)
            gi = _col(header, args.gene_col, "eQTL")
            vi = _col(header, args.variant_col, "eQTL")
            pi = header.index(args.p_col) if args.p_col in header else None
            if pi is not None:
                for line in fh:
                    f = _split(line, d)
                    if not f or f == [""]:
                        continue
                    key = (f[gi], f[vi])
                    if key[0] in gene_snps:
                        try:
                            pval[key] = float(f[pi])
                        except ValueError:
                            pass

    # ---- genotype: keep the cis universe, round to 0/1/2, standardise ------------------
    # genotype columns: variant_id, chromosome, position, ref, alt, then one per sample.
    with _open(args.genotype) as fh:
        geno_header = _split(fh.readline(), d)
        if geno_header[:5] != [args.variant_col, "chromosome", "position", "ref", "alt"]:
            raise SystemExit(
                "The genotype matrix must start with the columns "
                f"'{args.variant_col}', 'chromosome', 'position', 'ref', 'alt', "
                f"but it starts with {geno_header[:5]}.")
        samples = geno_header[5:]
        hardcall: dict[str, np.ndarray] = {}     # variant -> 0/1/2 vector over samples
        std_vec: dict[str, np.ndarray] = {}      # variant -> standardised vector (or absent)
        annot: dict[str, list[str]] = {}         # variant -> [chromosome, position, ref, alt]
        for line in fh:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            v = f[0]
            if v not in universe or v in hardcall:
                continue
            annot[v] = f[1:5]
            dosage = np.array(f[5:], dtype=np.float64)
            calls = np.clip(np.floor(dosage + 0.5), 0, 2).astype(np.uint8)
            hardcall[v] = calls
            sv = standardise(calls)
            if sv is not None:
                std_vec[v] = sv
    missing = universe - set(hardcall)
    if missing:
        raise SystemExit(f"{len(missing)} cis SNP(s) are absent from the genotype matrix, "
                         f"e.g. {sorted(missing)[:3]}.")

    # ---- prune each gene ---------------------------------------------------------------
    os.makedirs(args.outdir, exist_ok=True)
    pruned_path = os.path.join(args.outdir, f"{args.tissue}.cis_snps_pruned.tsv.gz")
    summ_path = os.path.join(args.outdir, f"{args.tissue}.ld_prune_summary.tsv")
    kept_universe: list[str] = []
    kept_seen: set[str] = set()
    n_in_total = n_kept_total = 0
    with _open(pruned_path, "wt") as fout, _open(summ_path, "wt") as fsum:
        fout.write(f"{args.gene_col}\t{args.variant_col}\n")
        fsum.write("gene_id\tn_cis_in\tn_cis_kept\tn_dropped_ld\tn_dropped_constant\n")
        for gene in sorted(gene_snps):
            snps = gene_snps[gene]
            if pval:
                order = sorted((s for s in snps if (gene, s) in pval),
                               key=lambda s: pval[(gene, s)])
                order += [s for s in snps if (gene, s) not in pval]   # unranked last, input order
            else:
                order = list(snps)
            kept, n_const, n_ld = greedy_clump(order, std_vec, args.ld_r2_threshold)
            for v in kept:
                fout.write(f"{gene}\t{v}\n")
                if v not in kept_seen:
                    kept_seen.add(v)
                    kept_universe.append(v)
            fsum.write(f"{gene}\t{len(snps)}\t{len(kept)}\t{n_ld}\t{n_const}\n")
            n_in_total += len(snps)
            n_kept_total += len(kept)

    # ---- hard-call genotype for the kept SNP universe ----------------------------------
    geno_path = os.path.join(args.outdir, f"{args.tissue}.genotype_012.tsv.gz")
    with _open(geno_path, "wt") as fout:
        fout.write(d.join([args.variant_col, "chromosome", "position", "ref", "alt"]
                          + samples) + "\n")
        for v in kept_universe:
            fout.write(d.join([v] + annot[v] + [str(int(x)) for x in hardcall[v]]) + "\n")

    # ---- QC report ---------------------------------------------------------------------
    qc_path = os.path.join(args.outdir, f"{args.tissue}.ld_prune_qc.tsv")
    qc = [
        ("ld_r2_threshold", args.ld_r2_threshold),
        ("genes", len(gene_snps)),
        ("samples", len(samples)),
        ("cis_snps_in_universe", len(universe)),
        ("cis_pairs_in", n_in_total),
        ("cis_pairs_kept", n_kept_total),
        ("kept_snp_universe", len(kept_universe)),
        ("ordered_by", "pvalue" if pval else "input_order"),
    ]
    with _open(qc_path, "wt") as fh:
        fh.write("metric\tvalue\n")
        for k, v in qc:
            fh.write(f"{k}\t{v}\n")
    sys.stderr.write(f"[ld_prune] tissue={args.tissue} genes={len(gene_snps)} "
                     f"pairs_in={n_in_total} pairs_kept={n_kept_total} "
                     f"kept_universe={len(kept_universe)} r2>{args.ld_r2_threshold}\n")


if __name__ == "__main__":
    main()
