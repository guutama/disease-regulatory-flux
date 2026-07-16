#!/usr/bin/env python3
"""
Select one genetic instrument per regulator gene and write the three matrices the GRN
reconstruction step (findr) consumes for one tissue.

A gene can act as a regulator only if it has a cis-eQTL. To orient a regulator->target
edge, findr uses a single genetic instrument for the regulator: its strongest cis-eQTL
(the "lead" SNP). This step reads one tissue's harmonised trio (the output of the
harmonisation step) and produces:

  --out-de   the instrument list: one row per regulator gene, giving its lead cis-eQTL as
             a (variant_id, gene_id) pair. The variant comes first and the gene second, the
             order findr expects (the SNP is the instrument, the gene is the regulator).

  --out-dg   the genotype of those instruments, with samples as rows and one column per
             instrument SNP (the dosage matrix findr reads, transposed from the
             variant-by-sample harmonised genotype).

  --out-dx   the expression of every gene, with samples as rows and one column per gene
             (transposed from the gene-by-sample harmonised expression). All genes are
             kept: every gene is a candidate target in the network, while only genes that
             appear in the instrument list act as regulators.

The lead cis-eQTL of a gene is the row with the smallest p-value. Ties are broken, in
order, by the largest absolute effect size and then the first variant id, so the choice is
deterministic and reproducible.

The expression (dX) and genotype (dG) matrices are written with their sample rows in the
same order, taken from the expression matrix, so the two line up row by row. The run stops
with an error if the expression and genotype matrices do not share the same samples, or if
a selected instrument is absent from the genotype matrix.

Inputs are the harmonised tables (tab-separated by default, optionally gzipped): the
expression matrix has the gene id in its first column and one column per sample; the
genotype matrix has the variant id in its first column and one column per sample; the eQTL
table has columns for the gene id, the variant id, the effect size, its standard error and
the p-value. The three output matrices are comma-separated, the format findr reads.

Column names and the input delimiter are options, so a new dataset is supported by adjusting
configuration rather than editing code.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from typing import TextIO


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


def _hardcall(dosage: str) -> str:
    """Round an allele dosage in [0, 2] to a 0/1/2 hard call. findr's causal tests need
    the instrument genotype as integer categories, not continuous (imputed) dosages."""
    return str(min(2, max(0, int(float(dosage) + 0.5))))


# ------------------------------------------------------------------- lead instrument choice
def select_leads(eqtl_rows, gene_i: int, var_i: int, beta_i: int, p_i: int):
    """Pick each gene's lead cis-eQTL from raw eQTL rows.

    eqtl_rows is an iterable of split string rows; gene_i/var_i/beta_i/p_i are the column
    indices of the gene id, variant id, effect size and p-value. The lead is the smallest
    p-value, ties broken by the largest |beta| and then the first variant id. Returns
    (leads, n_ties_broken) where leads maps gene_id -> variant_id, in first-seen gene order.
    """
    best: dict[str, tuple[float, float, str]] = {}   # gene -> (pvalue, -abs_beta, variant)
    leads: dict[str, str] = {}
    n_ties = 0
    for f in eqtl_rows:
        gene, variant = f[gene_i], f[var_i]
        pval = float(f[p_i])
        abs_beta = abs(float(f[beta_i]))
        key = (pval, -abs_beta, variant)
        if gene not in best:
            best[gene] = key
            leads[gene] = variant
            continue
        if pval == best[gene][0]:
            n_ties += 1
        if key < best[gene]:
            best[gene] = key
            leads[gene] = variant
    return leads, n_ties


# ----------------------------------------------------------------------------------- main
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # inputs (the harmonised trio)
    p.add_argument("--expression", required=True,
                   help="input: harmonised gene-by-sample expression matrix")
    p.add_argument("--genotype", required=True,
                   help="input: harmonised variant-by-sample dosage matrix")
    p.add_argument("--eqtl", required=True,
                   help="input: harmonised cis-eQTL table, one row per gene-variant pair")
    # outputs (the findr inputs)
    p.add_argument("--out-dx", required=True,
                   help="output: expression matrix, samples x genes (comma-separated)")
    p.add_argument("--out-dg", required=True,
                   help="output: instrument genotype, samples x SNPs (comma-separated)")
    p.add_argument("--out-de", required=True,
                   help="output: instrument list, variant_id,gene_id pairs (comma-separated)")
    p.add_argument("--qc", required=True,
                   help="output: TSV report of samples, genes and instruments selected")
    # format options (data-agnostic; wired from config)
    p.add_argument("--delimiter", default="\t",
                   help="field delimiter of the input tables (default: tab)")
    p.add_argument("--gene-col", default="gene_id",
                   help="gene-id column in the expression matrix and eQTL table "
                        "(default: gene_id)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column in the genotype matrix and eQTL table "
                        "(default: variant_id)")
    p.add_argument("--beta-col", default="beta",
                   help="eQTL effect-size column, used to break p-value ties (default: beta)")
    p.add_argument("--p-col", default="pvalue",
                   help="eQTL p-value column, used to pick the lead SNP (default: pvalue)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    # ---- eQTL: pick each gene's lead cis-eQTL (the instrument) -------------------------
    with _open(args.eqtl) as fh:
        header = _split(fh.readline(), d)
        idx = {name: i for i, name in enumerate(header)}
        for col in (args.gene_col, args.variant_col, args.beta_col, args.p_col):
            if col not in idx:
                raise SystemExit(f"The eQTL file is missing the required column '{col}'.")
        gi, vi, bi, pi = (idx[args.gene_col], idx[args.variant_col],
                          idx[args.beta_col], idx[args.p_col])
        rows = (f for f in (_split(line, d) for line in fh) if f and f != [""])
        leads, n_ties = select_leads(rows, gi, vi, bi, pi)
    if not leads:
        raise SystemExit("The eQTL table has no rows, so no instrument can be selected.")
    instrument_variants = set(leads.values())

    # ---- dG: genotype of the instruments, transposed to samples x SNPs -----------------
    # genotype columns: variant_id, chromosome, position, ref, alt, then one per sample.
    with _open(args.genotype) as fh:
        geno_header = _split(fh.readline(), d)
        if geno_header[:5] != [args.variant_col, "chromosome", "position", "ref", "alt"]:
            raise SystemExit(
                "The genotype matrix must start with the columns "
                f"'{args.variant_col}', 'chromosome', 'position', 'ref', 'alt', "
                f"but it starts with {geno_header[:5]}.")
        geno_samples = geno_header[5:]
        kept_variants: list[str] = []
        kept_dosages: list[list[str]] = []   # one inner list per instrument SNP, over geno_samples
        for line in fh:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            if f[0] in instrument_variants:
                kept_variants.append(f[0])
                kept_dosages.append(f[5:])
    missing = instrument_variants - set(kept_variants)
    if missing:
        raise SystemExit(f"{len(missing)} selected instrument(s) are absent from the genotype "
                         f"matrix, e.g. {sorted(missing)[:3]}. The genotype and eQTL files must "
                         "share these variants (they should after harmonisation).")

    # ---- dX: expression of all genes, transposed to samples x genes --------------------
    with _open(args.expression) as fh:
        expr_header = _split(fh.readline(), d)
        if expr_header[0] != args.gene_col:
            raise SystemExit(f"The first column of the expression matrix must be "
                             f"'{args.gene_col}', but it is '{expr_header[0]}'.")
        expr_samples = expr_header[1:]
        genes: list[str] = []
        expr_values: list[list[str]] = []    # one inner list per gene, over expr_samples
        for line in fh:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            genes.append(f[0])
            expr_values.append(f[1:])

    # ---- align genotype samples to the expression sample order -------------------------
    if set(expr_samples) != set(geno_samples):
        raise SystemExit("The expression and genotype matrices do not have the same samples. "
                         "They should after harmonisation; check the harmonised inputs.")
    geno_pos = {s: i for i, s in enumerate(geno_samples)}
    order = [geno_pos[s] for s in expr_samples]          # reorder genotype to expression order

    # ---- write the three findr matrices (comma-separated) ---------------------------
    with _open(args.out_dx, "wt") as fout:
        fout.write(",".join(genes) + "\n")
        for s in range(len(expr_samples)):
            fout.write(",".join(col[s] for col in expr_values) + "\n")

    with _open(args.out_dg, "wt") as fout:
        fout.write(",".join(kept_variants) + "\n")
        for s in order:
            fout.write(",".join(_hardcall(col[s]) for col in kept_dosages) + "\n")

    with _open(args.out_de, "wt") as fout:
        fout.write("variant_id,gene_id\n")
        for gene, variant in leads.items():
            fout.write(f"{variant},{gene}\n")

    # ---- QC report ---------------------------------------------------------------------
    qc = [
        ("samples", len(expr_samples)),
        ("genes_dX", len(genes)),
        ("regulators_dE", len(leads)),
        ("instruments_dG", len(kept_variants)),
        ("pvalue_ties_broken", n_ties),
    ]
    with _open(args.qc, "wt") as fh:
        fh.write("metric\tvalue\n")
        for k, v in qc:
            fh.write(f"{k}\t{v}\n")
    sys.stderr.write(f"[select_instruments] samples={len(expr_samples)} genes={len(genes)} "
                     f"regulators={len(leads)} instruments={len(kept_variants)} "
                     f"ties_broken={n_ties}\n")


if __name__ == "__main__":
    main()
