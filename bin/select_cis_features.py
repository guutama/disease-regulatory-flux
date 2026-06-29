#!/usr/bin/env python3
"""
Build each gene's set of cis-eQTL SNPs for one tissue --- the cis channel used to predict
the gene's expression.

The input is the harmonised cis-eQTL table, which already lists only cis associations (the
eQTL analysis restricted variants to each gene's cis window), so this step does not re-window
anything. It collects, per gene, the variants associated with it, dropping duplicate
gene-variant pairs. Optionally it keeps only the associations at or below a p-value
threshold; if the table has no p-value column the threshold is ignored and every association
is trusted as given.

Input:
  --eqtl   the harmonised cis-eQTL table, one row per gene-variant pair, with a gene-id and a
           variant-id column and, optionally, a p-value column.

Output (in --outdir), for tissue <tissue>:
  <tissue>.cis_features.tsv.gz       one row per kept pair: gene_id, variant_id;
  <tissue>.cis_features_summary.tsv  per gene: number of cis SNPs kept;
  <tissue>.cis_features_qc.tsv       pairs read, kept and dropped, and whether the p-value
                                     filter was applied.

Column names and the field delimiter are options, so a new dataset is supported by adjusting
configuration rather than editing code.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict
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


# ----------------------------------------------------------------------------------- main
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eqtl", required=True,
                   help="input: harmonised cis-eQTL table (gene-variant pairs)")
    p.add_argument("--tissue", required=True, help="tissue label, used to name the outputs")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--max-pvalue", type=float, default=None,
                   help="if set and the eQTL table has a p-value column, keep only pairs "
                        "with p-value <= this; ignored when no p-value column is present")
    # format options
    p.add_argument("--delimiter", default="\t",
                   help="field delimiter of the eQTL table (default: tab)")
    p.add_argument("--gene-col", default="gene_id",
                   help="gene-id column in the eQTL table (default: gene_id)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column in the eQTL table (default: variant_id)")
    p.add_argument("--p-col", default="pvalue",
                   help="p-value column in the eQTL table, if present (default: pvalue)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    with _open(args.eqtl) as fh:
        header = _split(fh.readline(), d)
        idx = {name: i for i, name in enumerate(header)}
        for col in (args.gene_col, args.variant_col):
            if col not in idx:
                raise SystemExit(f"The eQTL file is missing the required column '{col}'.")
        gi, vi = idx[args.gene_col], idx[args.variant_col]
        pi = idx.get(args.p_col)            # None if the table has no p-value column
        apply_p = args.max_pvalue is not None and pi is not None

        cis: dict[str, list[str]] = defaultdict(list)
        seen: set[tuple[str, str]] = set()
        n_in = n_dropped_p = 0
        for line in fh:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            n_in += 1
            if apply_p and float(f[pi]) > args.max_pvalue:
                n_dropped_p += 1
                continue
            gene, variant = f[gi], f[vi]
            if (gene, variant) not in seen:        # keep first-seen order, drop duplicates
                seen.add((gene, variant))
                cis[gene].append(variant)

    if not cis:
        raise SystemExit("No cis-eQTL pair remains, so no cis feature set can be built.")

    os.makedirs(args.outdir, exist_ok=True)
    feat_path = os.path.join(args.outdir, f"{args.tissue}.cis_features.tsv.gz")
    summ_path = os.path.join(args.outdir, f"{args.tissue}.cis_features_summary.tsv")
    n_kept = 0
    with _open(feat_path, "wt") as fout, _open(summ_path, "wt") as fsum:
        fout.write("gene_id\tvariant_id\n")
        fsum.write("gene_id\tn_cis_snps\n")
        for gene in sorted(cis):
            variants = cis[gene]
            for variant in variants:
                fout.write(f"{gene}\t{variant}\n")
                n_kept += 1
            fsum.write(f"{gene}\t{len(variants)}\n")

    qc_path = os.path.join(args.outdir, f"{args.tissue}.cis_features_qc.tsv")
    qc = [
        ("genes", len(cis)),
        ("pairs_in", n_in),
        ("pairs_kept", n_kept),
        ("dropped_pvalue", n_dropped_p),
        ("dropped_duplicate", n_in - n_dropped_p - n_kept),
        ("pvalue_filter_applied", "yes" if apply_p else "no"),
    ]
    with _open(qc_path, "wt") as fh:
        fh.write("metric\tvalue\n")
        for k, v in qc:
            fh.write(f"{k}\t{v}\n")
    sys.stderr.write(f"[cis_features] tissue={args.tissue} genes={len(cis)} "
                     f"pairs_kept={n_kept} dropped_pvalue={n_dropped_p} "
                     f"pvalue_filter={'yes' if apply_p else 'no'}\n")


if __name__ == "__main__":
    main()
