#!/usr/bin/env python3
"""
Shared gene-label utility: map a versioned-ENSG ``gene_id`` to its HGNC symbol.

Single source of truth for presenting results by human-readable gene name
(e.g. ``ENSG00000118526.6`` -> ``TCF21``) while keeping ``gene_id`` as the
unambiguous key. Used two ways:

  * as a library (``import gene_labels``) by the plotting / analysis scripts:
        m = gene_labels.load_label_map(annot_path)
        df = gene_labels.add_symbol_column(df, mapping=m)        # adds 'gene'
        label = gene_labels.resolve("ENSG00000118526.6", m)      # -> 'TCF21'

  * as a CLI to batch-relabel an existing result TSV in place / to a new file:
        gene_labels.py --in scores.tsv --out scores.labeled.tsv \\
            --gene-annot gencode.v19.genes.tsv

Resolution mirrors bin/harmonise_expression.py exactly: try the versioned
ENSG, then the version-stripped base ENSG, else fall back to the base ENSG
itself (never empty). Because several ENSG can share a symbol, the symbol is
NOT unique -- always keep ``gene_id`` alongside it as the disambiguator.

The annotation file is the project's GENCODE v19 (GRCh37) table with columns
``gene_id`` (versioned ENSG) and ``gene_name`` (HGNC symbol); the same file
Stage 2 uses. Path is data-agnostic: pass --gene-annot, or set $GENE_ANNOT.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Optional


def strip_ensg_version(gene_id: str) -> str:
    """ENSG00000254681.2 -> ENSG00000254681. Unchanged if there is no dot."""
    dot = gene_id.find(".")
    return gene_id[:dot] if dot >= 0 else gene_id


def load_label_map(path: str, id_col: str = "gene_id",
                   symbol_col: str = "gene_name") -> dict[str, str]:
    """Build a gene_id -> HGNC symbol dict keyed by BOTH the versioned ENSG and
    its base (version-stripped) form, so either lookup form resolves directly."""
    m: dict[str, str] = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_id, i_sym = header.index(id_col), header.index(symbol_col)
        except ValueError as e:
            sys.exit(f"ERROR: gene-annotation column missing: {e}; header={header}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_id, i_sym):
                continue
            gid, name = f[i_id], f[i_sym]
            if gid and name:
                m[gid] = name
                m.setdefault(strip_ensg_version(gid), name)
    if not m:
        sys.exit(f"ERROR: gene-annotation map is empty: {path}")
    return m


def resolve(gene_id: str, mapping: dict[str, str]) -> str:
    """versioned ENSG -> HGNC symbol, with version-strip fallback, then base ENSG."""
    s = mapping.get(gene_id)
    if s:
        return s
    base = strip_ensg_version(gene_id)
    return mapping.get(base, base)


def add_symbol_column(df, mapping: dict[str, str], id_col: str = "gene_id",
                      out_col: str = "gene", position: int = 0):
    """Return ``df`` with an HGNC-symbol column (default name 'gene') inserted at
    ``position`` (default: first). ``id_col`` is preserved as the disambiguator.
    Idempotent: an existing ``out_col`` is overwritten in place."""
    sym = df[id_col].map(lambda g: resolve(str(g), mapping))
    if out_col in df.columns:
        df[out_col] = sym
        return df
    df = df.copy()
    df.insert(position, out_col, sym)
    return df


def default_annot(explicit: Optional[str]) -> str:
    path = explicit or os.environ.get("GENE_ANNOT")
    if not path:
        sys.exit("ERROR: no gene annotation given (use --gene-annot or $GENE_ANNOT)")
    if not os.path.exists(path):
        sys.exit(f"ERROR: gene annotation not found: {path}")
    return path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in", dest="infile", required=True, help="input TSV (keyed on --id-col)")
    p.add_argument("--out", required=True, help="output TSV with the symbol column added")
    p.add_argument("--gene-annot", default=None, help="GENCODE TSV (or set $GENE_ANNOT)")
    p.add_argument("--id-col", default="gene_id")
    p.add_argument("--out-col", default="gene")
    p.add_argument("--gene-annot-id-col", default="gene_id")
    p.add_argument("--gene-annot-symbol-col", default="gene_name")
    p.add_argument("--sep", default="\t")
    a = p.parse_args()

    import pandas as pd
    mapping = load_label_map(default_annot(a.gene_annot),
                             a.gene_annot_id_col, a.gene_annot_symbol_col)
    df = pd.read_csv(a.infile, sep=a.sep)
    if a.id_col not in df.columns:
        sys.exit(f"ERROR: --id-col '{a.id_col}' not in {a.infile} (cols: {list(df.columns)})")
    df = add_symbol_column(df, mapping, id_col=a.id_col, out_col=a.out_col)
    df.to_csv(a.out, sep=a.sep, index=False)
    n_unmapped = int((df[a.out_col] == df[a.id_col].map(strip_ensg_version)).sum())
    print(f"wrote {a.out}: {len(df):,} rows, {n_unmapped:,} unmapped (symbol = base ENSG)")


if __name__ == "__main__":
    main()
