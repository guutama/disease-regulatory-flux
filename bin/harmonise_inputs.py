#!/usr/bin/env python3
"""
Align the three reference inputs of the disease-regulatory flux pipeline onto one shared
set of samples.

The method models every expressed gene as a potential target in the regulatory network,
so all genes in the expression matrix are kept. Only genes that have a cis-eQTL can act as
regulators, and those associations come from the eQTL map. The genotype is used only for
the variants that appear in the eQTL map, so it is reduced to those variants. Samples are
therefore the one thing that has to be matched across the files.

Provide one tissue's three reference files:

  --expression   Gene-by-sample matrix of analysis-ready expression (already normalised
                 and covariate-adjusted). The first column holds the gene id; the
                 remaining columns are the samples (one column per sample), which should
                 be the same samples as in the genotype matrix.

  --genotype     Variant-by-sample matrix of alternate-allele dosages in [0, 2]. The
                 leading column(s) identify the variant (see "Variant identifiers"); the
                 remaining columns are the samples (one column per sample).

  --eqtl         Cis-eQTL table with one row per gene-variant pair, with columns for the
                 gene id, the variant identifier, the effect size (beta), its standard
                 error (se) and the p-value.

The three files are linked by three keys: the sample columns connect the expression and
genotype matrices; the gene id connects the expression matrix and the eQTL table; the
variant identifier connects the genotype matrix and the eQTL table.

Samples
  The expression and genotype matrices must share samples. Both are reduced to the samples
  they have in common and written in one consistent order. The run stops with an error if
  they share no samples.

Variant identifiers
  The genotype and eQTL files must identify variants the same way. If both files contain
  an rsID column, variants are matched by rsID. Otherwise both files must provide the
  chromosome, position, ref and alt columns, and variants are matched by
  "chr:pos:ref:alt". The run stops with an error if the two files share neither form.

Outputs
  - expression: every gene, with the columns reduced and reordered to the shared samples;
  - genotype:   reduced to the variants that appear in the eQTL map, with the same shared
                samples in the same order;
  - eQTL:       the cis pairs whose variant is present in the genotype and whose gene is
                present in the expression matrix;
  - a QC report of how many samples, genes and variants were kept and dropped.
  Expression values and dosages are copied unchanged; only rows and columns are selected.
  The run stops with an error if no eQTL pair survives the matching.

Column names and the field delimiter are options, so a new dataset is supported by
adjusting configuration rather than editing code.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from typing import Callable, TextIO


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


# -------------------------------------------------------------------- variant key handling
def resolve_key_type(geno_header: list[str], eqtl_header: list[str], args) -> str:
    """Choose the variant key shared by genotype and eQTL: 'rsid' if both carry the rsID
    column, else 'composite' if both carry chrom/pos/ref/alt, else raise SystemExit."""
    g, e = set(geno_header), set(eqtl_header)
    if args.rsid_col in g and args.rsid_col in e:
        return "rsid"
    composite = {args.chrom_col, args.pos_col, args.ref_col, args.alt_col}
    if composite <= g and composite <= e:
        return "composite"
    raise SystemExit(
        "The genotype and eQTL files do not identify variants the same way.\n"
        f"  Provide an '{args.rsid_col}' column in both files, or provide the "
        f"{args.chrom_col}/{args.pos_col}/{args.ref_col}/{args.alt_col} columns in both.\n"
        f"  genotype columns: {geno_header}\n  eQTL columns:     {eqtl_header}"
    )


def variant_id_getter(header: list[str], key_type: str, args) -> Callable[[list[str]], str]:
    """Return a function row->variant_id for a file with this header under key_type."""
    idx = {name: i for i, name in enumerate(header)}
    if key_type == "rsid":
        j = idx[args.rsid_col]
        return lambda row: row[j]
    c, p, r, a = idx[args.chrom_col], idx[args.pos_col], idx[args.ref_col], idx[args.alt_col]
    return lambda row: f"{row[c]}:{row[p]}:{row[r]}:{row[a]}"


def id_columns(key_type: str, args) -> set[str]:
    """The non-sample identifier columns of the genotype matrix: the variant key plus
    the chrom/pos/ref/alt allele annotation that is always carried alongside it."""
    alleles = {args.chrom_col, args.pos_col, args.ref_col, args.alt_col}
    if key_type == "rsid":
        return {args.rsid_col} | alleles
    return alleles


# ------------------------------------------------------------------------------ readers
def read_header(path: str, delim: str) -> list[str]:
    with _open(path) as fh:
        return _split(fh.readline(), delim)


def read_gene_set(path: str, delim: str) -> set[str]:
    """Genes present in the expression matrix (the first column, excluding the header)."""
    genes: set[str] = set()
    with _open(path) as fh:
        fh.readline()
        for line in fh:
            f = _split(line, delim)
            if f and f[0] != "":
                genes.add(f[0])
    return genes


def read_eqtl(path: str, args, key_type: str):
    """Read the (small) eQTL table fully. Returns (rows, genes, variants) where rows is a
    list of (gene_id, variant_id, beta, se, pvalue) raw-string tuples."""
    rows: list[tuple[str, str, str, str, str]] = []
    genes: set[str] = set()
    variants: set[str] = set()
    with _open(path) as fh:
        header = _split(fh.readline(), args.delimiter)
        idx = {name: i for i, name in enumerate(header)}
        for col in (args.gene_col, args.beta_col, args.se_col, args.p_col):
            if col not in idx:
                raise SystemExit(f"The eQTL file is missing the required column '{col}'.")
        get_vid = variant_id_getter(header, key_type, args)
        gi, bi, si, pi = idx[args.gene_col], idx[args.beta_col], idx[args.se_col], idx[args.p_col]
        for line in fh:
            f = _split(line, args.delimiter)
            if not f or f == [""]:
                continue
            gene, vid = f[gi], get_vid(f)
            rows.append((gene, vid, f[bi], f[si], f[pi]))
            genes.add(gene)
            variants.add(vid)
    return rows, genes, variants


# ----------------------------------------------------------------------------------- main
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # inputs
    p.add_argument("--expression", required=True,
                   help="input: gene-by-sample expression matrix (gene id in column 1)")
    p.add_argument("--genotype", required=True,
                   help="input: variant-by-sample dosage matrix (dosages in [0, 2])")
    p.add_argument("--eqtl", required=True,
                   help="input: cis-eQTL table, one row per gene-variant pair")
    # outputs
    p.add_argument("--out-expression", required=True,
                   help="output: expression matrix aligned to the shared samples")
    p.add_argument("--out-genotype", required=True,
                   help="output: genotype reduced to the eQTL variants and shared samples")
    p.add_argument("--out-eqtl", required=True,
                   help="output: eQTL table restricted to the matched genes and variants")
    p.add_argument("--qc", required=True,
                   help="output: TSV report of samples/genes/variants kept and dropped")
    # format options (data-agnostic; wired from config)
    p.add_argument("--delimiter", default="\t",
                   help="field delimiter of all input tables (default: tab)")
    p.add_argument("--gene-col", default="gene_id",
                   help="gene-id column in the expression matrix and eQTL table "
                        "(default: gene_id)")
    p.add_argument("--rsid-col", default="rsID",
                   help="rsID variant column, if present in genotype and eQTL "
                        "(default: rsID)")
    p.add_argument("--chrom-col", default="chromosome",
                   help="chromosome column, used for the chr:pos:ref:alt variant key "
                        "(default: chromosome)")
    p.add_argument("--pos-col", default="position",
                   help="position column for the chr:pos:ref:alt key (default: position)")
    p.add_argument("--ref-col", default="ref",
                   help="reference-allele column for the chr:pos:ref:alt key (default: ref)")
    p.add_argument("--alt-col", default="alt",
                   help="alternate-allele column for the chr:pos:ref:alt key (default: alt)")
    p.add_argument("--beta-col", default="beta",
                   help="eQTL effect-size column (default: beta)")
    p.add_argument("--se-col", default="se",
                   help="eQTL standard-error column (default: se)")
    p.add_argument("--p-col", default="pvalue",
                   help="eQTL p-value column (default: pvalue)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    # ---- headers, key type, eQTL table -------------------------------------------------
    expr_header = read_header(args.expression, d)
    if expr_header[0] != args.gene_col:
        raise SystemExit(f"The first column of the expression matrix must be '{args.gene_col}', "
                         f"but it is '{expr_header[0]}'.")
    expr_samples = expr_header[1:]

    geno_header = read_header(args.genotype, d)
    eqtl_header = read_header(args.eqtl, d)
    key_type = resolve_key_type(geno_header, eqtl_header, args)

    # chrom/pos/ref/alt must be present so each variant can later be aligned to a GWAS
    # by its alleles, whatever variant key is used above. They are carried through to
    # the genotype output and on to the downstream stages.
    allele_cols = (args.chrom_col, args.pos_col, args.ref_col, args.alt_col)
    missing_allele = [c for c in allele_cols if c not in geno_header]
    if missing_allele:
        raise SystemExit(
            "The genotype matrix must carry the allele-annotation columns "
            f"{list(allele_cols)} so variants can be aligned to a GWAS downstream.\n"
            f"  missing: {missing_allele}\n  genotype columns: {geno_header}")

    eqtl_rows, eqtl_genes, eqtl_variants = read_eqtl(args.eqtl, args, key_type)
    expr_genes = read_gene_set(args.expression, d)

    geno_id_cols = id_columns(key_type, args)
    geno_samples = [c for c in geno_header if c not in geno_id_cols]

    # ---- align on samples (order follows the expression matrix) ------------------------
    geno_sample_set = set(geno_samples)
    keep_samples = [s for s in expr_samples if s in geno_sample_set]
    if not keep_samples:
        raise SystemExit("The expression and genotype matrices share no samples. "
                         "Check that their sample column names match.")

    # ---- stream genotype: keep variants present in eQTL, reorder to keep_samples -------
    get_vid = variant_id_getter(geno_header, key_type, args)
    geno_idx = {c: i for i, c in enumerate(geno_header)}
    keep_col_idx = [geno_idx[s] for s in keep_samples]
    ci, pi, ri, ai = (geno_idx[args.chrom_col], geno_idx[args.pos_col],
                      geno_idx[args.ref_col], geno_idx[args.alt_col])
    kept_variants: set[str] = set()
    n_geno_total = 0
    with _open(args.genotype) as fin, _open(args.out_genotype, "wt") as fout:
        fin.readline()
        fout.write(d.join(["variant_id", "chromosome", "position", "ref", "alt"]
                          + keep_samples) + "\n")
        for line in fin:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            n_geno_total += 1
            vid = get_vid(f)
            if vid not in eqtl_variants:
                continue
            kept_variants.add(vid)
            fout.write(d.join([vid, f[ci], f[pi], f[ri], f[ai]]
                              + [f[i] for i in keep_col_idx]) + "\n")
    if not kept_variants:
        raise SystemExit("The genotype matrix and the eQTL table share no variants. "
                         "Check that their variant identifiers match.")

    # ---- eQTL: cis pairs whose variant is in the genotype and whose gene is expressed --
    n_eqtl = 0
    regulators: set[str] = set()
    with _open(args.out_eqtl, "wt") as fout:
        fout.write(d.join(["gene_id", "variant_id", "beta", "se", "pvalue"]) + "\n")
        for gene, vid, beta, se, pval in eqtl_rows:
            if vid in kept_variants and gene in expr_genes:
                fout.write(d.join([gene, vid, beta, se, pval]) + "\n")
                regulators.add(gene)
                n_eqtl += 1
    if not n_eqtl:
        raise SystemExit("No cis-eQTL pair remains after matching to the genotype and "
                         "expression. Check that variant ids and gene ids match across files.")

    # ---- aligned expression: ALL genes, columns reduced/reordered to shared samples ----
    expr_idx = {c: i for i, c in enumerate(expr_header)}
    expr_keep_idx = [expr_idx[s] for s in keep_samples]
    n_expr = 0
    with _open(args.expression) as fin, _open(args.out_expression, "wt") as fout:
        fin.readline()
        fout.write(d.join([args.gene_col] + keep_samples) + "\n")
        for line in fin:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            fout.write(d.join([f[0]] + [f[i] for i in expr_keep_idx]) + "\n")
            n_expr += 1

    # ---- QC report ---------------------------------------------------------------------
    qc = [
        ("variant_key_type", key_type),
        ("samples_expression", len(expr_samples)),
        ("samples_genotype", len(geno_samples)),
        ("samples_kept", len(keep_samples)),
        ("genes_expression", len(expr_genes)),
        ("genes_written", n_expr),
        ("genes_eqtl", len(eqtl_genes)),
        ("regulators_kept", len(regulators)),
        ("variants_genotype", n_geno_total),
        ("variants_eqtl", len(eqtl_variants)),
        ("variants_kept", len(kept_variants)),
        ("eqtl_pairs_total", len(eqtl_rows)),
        ("eqtl_pairs_kept", n_eqtl),
    ]
    with _open(args.qc, "wt") as fh:
        fh.write("metric\tvalue\n")
        for k, v in qc:
            fh.write(f"{k}\t{v}\n")
    sys.stderr.write(f"[harmonise] samples={len(keep_samples)} genes={n_expr} "
                     f"regulators={len(regulators)} variants={len(kept_variants)} "
                     f"eqtl_pairs={n_eqtl} key={key_type}\n")


if __name__ == "__main__":
    main()
