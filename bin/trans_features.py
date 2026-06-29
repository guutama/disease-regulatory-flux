#!/usr/bin/env python3
"""
Derive each gene's trans-feature SNPs for one tissue by walking the regulatory network
upward from the gene to its regulators.

A gene's expression is predicted from two genetic channels: its own cis-eQTLs (the cis
channel) and the cis-eQTLs of the regulators upstream of it in the network (the trans
channel). This step builds the trans channel. Starting from each gene, it follows the
network edges backwards --- to the gene's direct regulators (hop 1), their regulators
(hop 2), and so on up to a chosen number of hops --- and collects, for every regulator
reached, that regulator's LD-pruned cis SNPs. Those SNPs are the gene's trans features.

Because every regulator in the network is anchored by its own cis-eQTL, every ancestor
reached on the walk has a cis-SNP set to contribute. A regulator that is reachable at more
than one hop (along paths of different length) is reported at each such hop, so a downstream
model can keep or collapse hops as it chooses.

Inputs:
  --grn          the tissue's reconstructed network (the GRN reconstruction step), with a
                 source-gene and a target-gene column; an edge source -> target means the
                 source regulates the target.
  --cis-pruned   the LD-pruned per-gene cis-SNP set (the LD-filtering step): gene_id,
                 variant_id. This gives each regulator its set of cis SNPs.

Output (in --outdir), for tissue <tissue>:
  <tissue>.trans_features.tsv.gz       one row per (gene, hop, regulator, variant):
                                       gene_id, hop, source_gene_id, variant_id;
  <tissue>.trans_features_summary.tsv  per gene: number of regulators and unique trans SNPs;
  <tissue>.trans_features_qc.tsv       tissue-level totals.

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


# -------------------------------------------------------------------- network walk upward
def ancestors_by_hop(parents: dict[str, list[str]], gene: str, max_hop: int):
    """Yield (hop, ancestor) for `gene`: at each hop 1..max_hop, the genes reachable from
    `gene` in exactly that many backward steps. A gene reachable at several hop counts (via
    paths of different length) is yielded at each."""
    current = {gene}
    for hop in range(1, max_hop + 1):
        nxt: set[str] = set()
        for node in current:
            nxt.update(parents.get(node, ()))
        if not nxt:
            return
        for anc in sorted(nxt):
            yield hop, anc
        current = nxt


# ----------------------------------------------------------------------------------- main
def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--grn", required=True,
                   help="input: reconstructed network edge table (source/target columns)")
    p.add_argument("--cis-pruned", required=True,
                   help="input: LD-pruned per-gene cis-SNP set (gene_id, variant_id)")
    p.add_argument("--tissue", required=True, help="tissue label, used to name the outputs")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--max-hop", type=int, default=2,
                   help="how many regulator hops to walk upward (default: 2)")
    # format options
    p.add_argument("--grn-delimiter", default=",",
                   help="field delimiter of the GRN edge table (default: comma)")
    p.add_argument("--cis-delimiter", default="\t",
                   help="field delimiter of the cis-SNP table (default: tab)")
    p.add_argument("--source-col", default="Source",
                   help="regulator (source-gene) column in the GRN table (default: Source)")
    p.add_argument("--target-col", default="Target",
                   help="target-gene column in the GRN table (default: Target)")
    p.add_argument("--gene-col", default="gene_id",
                   help="gene-id column in the cis-SNP table (default: gene_id)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column in the cis-SNP table (default: variant_id)")
    a = p.parse_args(argv)
    for attr in ("grn_delimiter", "cis_delimiter"):
        if getattr(a, attr) == "\\t":
            setattr(a, attr, "\t")
    if a.max_hop < 1:
        raise SystemExit("--max-hop must be at least 1.")
    return a


def main(argv=None) -> None:
    args = parse_args(argv)

    # ---- pruned cis-SNP set of every regulator gene ------------------------------------
    cis: dict[str, list[str]] = defaultdict(list)
    with _open(args.cis_pruned) as fh:
        header = _split(fh.readline(), args.cis_delimiter)
        gi = _col(header, args.gene_col, "cis-SNP")
        vi = _col(header, args.variant_col, "cis-SNP")
        for line in fh:
            f = _split(line, args.cis_delimiter)
            if not f or f == [""]:
                continue
            cis[f[gi]].append(f[vi])

    # ---- child -> parents map from the network edges -----------------------------------
    parents: dict[str, list[str]] = defaultdict(list)
    genes: set[str] = set()
    with _open(args.grn) as fh:
        header = _split(fh.readline(), args.grn_delimiter)
        si = _col(header, args.source_col, "GRN")
        ti = _col(header, args.target_col, "GRN")
        for line in fh:
            f = _split(line, args.grn_delimiter)
            if not f or f == [""]:
                continue
            source, target = f[si], f[ti]
            parents[target].append(source)
            genes.add(source)
            genes.add(target)

    # ---- walk each gene upward and emit its trans features -----------------------------
    os.makedirs(args.outdir, exist_ok=True)
    feat_path = os.path.join(args.outdir, f"{args.tissue}.trans_features.tsv.gz")
    summ_path = os.path.join(args.outdir, f"{args.tissue}.trans_features_summary.tsv")
    n_rows = 0
    n_genes_with_features = 0
    with _open(feat_path, "wt") as fout, _open(summ_path, "wt") as fsum:
        fout.write("gene_id\thop\tsource_gene_id\tvariant_id\n")
        fsum.write("gene_id\tn_regulators\tn_unique_snps\n")
        for gene in sorted(genes):
            regulators: set[str] = set()
            snps: set[str] = set()
            for hop, anc in ancestors_by_hop(parents, gene, args.max_hop):
                anc_cis = cis.get(anc)
                if not anc_cis:
                    continue
                regulators.add(anc)
                for variant in anc_cis:
                    fout.write(f"{gene}\t{hop}\t{anc}\t{variant}\n")
                    snps.add(variant)
                    n_rows += 1
            if regulators:
                n_genes_with_features += 1
                fsum.write(f"{gene}\t{len(regulators)}\t{len(snps)}\n")

    # ---- QC report ---------------------------------------------------------------------
    qc_path = os.path.join(args.outdir, f"{args.tissue}.trans_features_qc.tsv")
    qc = [
        ("max_hop", args.max_hop),
        ("genes_in_network", len(genes)),
        ("regulators_with_cis", sum(1 for g in genes if cis.get(g))),
        ("genes_with_trans_features", n_genes_with_features),
        ("trans_feature_rows", n_rows),
    ]
    with _open(qc_path, "wt") as fh:
        fh.write("metric\tvalue\n")
        for k, v in qc:
            fh.write(f"{k}\t{v}\n")
    sys.stderr.write(f"[trans_features] tissue={args.tissue} genes={len(genes)} "
                     f"with_features={n_genes_with_features} rows={n_rows} "
                     f"max_hop={args.max_hop}\n")


if __name__ == "__main__":
    main()
