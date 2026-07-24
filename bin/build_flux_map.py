#!/usr/bin/env python3
"""Step 10 -- disease-regulatory flux map: a path-preserving network of the regulators that
deliver trans signal to each disease gene.

The trans component of a gene's association (Step 9) is, by construction of the trans channel,
an exact sum of contributions from the gene's upstream GRN regulators: every trans SNP is a
cis-eQTL of one or more of the gene's ancestors. For a disease gene g and a regulator A that
feeds it, the flux A delivers is

    flux(A -> g) = (1/sigma_g) * sum_{s in E(A) cap T(g)} w_{g,s} * sd_s * z_gwas_s ,

where E(A) is A's cis-eQTL set and T(g) is g's trans-feature set; a SNP that instruments several
ancestors splits its term evenly among them, so sum_A flux(A -> g) = z_trans(g).

The map keeps the network paths explicit rather than collapsing them. Its nodes are the disease
genes (significant at p_adj < --fdr) together with their hop-1 GRN parents and hop-2
grandparents. Its edges are the true single-hop GRN edges along the paths down to a disease
gene -- recovered from the hop-1 trans features (source_gene -> gene):

  * parent -> disease_gene            (hop 1)
  * grandparent -> intermediate_parent (hop 2; the intermediate parent is kept, not collapsed)

Each edge carries the signed flux its source regulator delivers to the disease gene it feeds
(positive = drives the gene toward the risk-increasing direction). A grandparent that reaches a
disease gene through several parents splits its flux evenly across those edges, so the fluxes
toward a gene still sum to z_trans(g). Every parent and grandparent is kept even when it delivers
no flux (its `flux` is left blank). This is an algebraic decomposition, not a fitted model.

Inputs
  --association     association_<tissue>_<trait>.tsv (Step 9): gene_id, z_trans, sigma_g, p_adj
  --weights         <tissue>.expr_model_weights.tsv.gz (Step 7): the trans-channel weights
  --trans-features  <tissue>.trans_features.tsv.gz (Step 6): gene_id, hop, source_gene_id, variant_id
  --genotype        <tissue>.genotype_012.tsv.gz (Step 5): supplies sd_s
  --gwas            <trait>.gwas.tsv.gz (Step 8): ALT-aligned z
  --tissue --trait  labels    --outdir    [--fdr 0.05]

Output: <outdir>/flux_<tissue>_<trait>.csv, one row per regulator->disease-gene attribution
  regulator  regulator_gene_id  target  target_gene_id  tissue  hop  flux
Each contributing ancestor is attributed DIRECTLY to the disease gene it feeds: hop=1 if it is a
direct GRN parent, hop=2 if it is a grandparent (its hop-2 shares summed back to the disease gene).
An ancestor that is both a parent and a grandparent of the same gene is kept once at hop=1, so the
fluxes into a target still sum to z_trans(g). regulator/target are HGNC symbols; the *_gene_id
columns carry the (versioned Ensembl) identifiers. This is the manuscript Supplementary Table S4
column layout.
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import gene_labels as gl  # noqa: E402

DEFAULT_ANNOT = ("/cluster/projects/nn1015k/datasets/references/gencode/GRCh37/"
                 "gencode.v19.genes.tsv")


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


def load_association(path, delim, gene_col, fdr):
    """gene_id -> (sigma_g, z_trans) for genes significant at p_adj < fdr."""
    sig: dict[str, tuple[float, float]] = {}
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        gi = _col(header, gene_col, "association")
        zi = _col(header, "z_trans", "association")
        si = _col(header, "sigma_g", "association")
        pi = _col(header, "p_adj", "association")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            if float(f[pi]) < fdr:
                sig[f[gi]] = (float(f[si]), float(f[zi]))
    return sig


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


def load_trans_weights(path, delim, gene_col, variant_col, selected):
    """gene_id -> {variant_id: weight} for the trans channel of each gene's SELECTED config.

    The weights file holds every fitted configuration; we keep only the trans rows of the
    winning config, so the flux uses the same predictor the association did."""
    w: dict[str, dict[str, float]] = defaultdict(dict)
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
            if f[chi] == "trans" and selected.get(f[gi]) == f[ci]:
                w[f[gi]][f[vi]] = float(f[wi])
    return w


def load_trans_features(path, delim, gene_col, variant_col):
    """Returns (per_gene, parents). per_gene[g] is a list of (hop, source_gene, variant);
    parents[g] is the set of g's hop-1 GRN parents (the single-hop GRN, source -> gene)."""
    per_gene: dict[str, list[tuple[str, str, str]]] = defaultdict(list)
    parents: dict[str, set] = defaultdict(set)
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        gi = _col(header, gene_col, "trans-features")
        hi = _col(header, "hop", "trans-features")
        si = _col(header, "source_gene_id", "trans-features")
        vi = _col(header, variant_col, "trans-features")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            per_gene[f[gi]].append((f[hi], f[si], f[vi]))
            if f[hi] == "1":
                parents[f[gi]].add(f[si])
    return per_gene, parents


def load_genotype_sd(path, delim, variant_col, wanted):
    """variant_id -> genotype standard deviation, for the wanted variants."""
    sd: dict[str, float] = {}
    with _open(path) as fh:
        header = _split(fh.readline(), delim)
        if header[:5] != [variant_col, "chromosome", "position", "ref", "alt"]:
            raise SystemExit(
                f"The genotype matrix must start with '{variant_col}', 'chromosome', "
                f"'position', 'ref', 'alt', but it starts with {header[:5]}.")
        for line in fh:
            f = _split(line, delim)
            if not f or f == [""]:
                continue
            if f[0] in wanted and f[0] not in sd:
                sd[f[0]] = float(np.asarray(f[5:], dtype=np.float64).std())
    return sd


def load_gwas_z(path, delim, variant_col):
    """variant_id -> ALT-aligned z."""
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


def gene_flux(rows, wmap, sd, gwas_z, sigma_g):
    """For one disease gene, return source_gene -> (flux, n_snp) over its trans features.

    A SNP shared by several ancestors splits its term evenly; an ancestor with no weighted,
    GWAS-matched SNP gets (0.0, 0)."""
    ancestors_per_snp = defaultdict(set)
    for _hop, A, s in rows:
        ancestors_per_snp[s].add(A)
    acc = defaultdict(lambda: [0.0, 0])
    seen = set()
    for _hop, A, s in rows:
        if (A, s) in seen:
            continue
        seen.add((A, s))
        w = wmap.get(s)
        z = gwas_z.get(s)
        ss = sd.get(s)
        if w is None or z is None or ss is None:
            continue
        acc[A][0] += w * ss * z / sigma_g / len(ancestors_per_snp[s])
        acc[A][1] += 1
    return acc


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--association", required=True, help="input: Step 9 association table")
    p.add_argument("--weights", required=True, help="input: expr_model_weights matrix")
    p.add_argument("--selected", required=True,
                   help="input: expr_model_selected table (each predictable gene's best config)")
    p.add_argument("--trans-features", required=True, help="input: trans_features table")
    p.add_argument("--genotype", required=True, help="input: genotype_012 hard-call matrix")
    p.add_argument("--gwas", required=True, help="input: ALT-aligned GWAS (Step 8)")
    p.add_argument("--tissue", required=True, help="tissue label")
    p.add_argument("--trait", required=True, help="trait label")
    p.add_argument("--outdir", required=True, help="output directory")
    p.add_argument("--fdr", type=float, default=0.05,
                   help="p_adj threshold for a gene to be a disease gene (default: 0.05)")
    p.add_argument("--delimiter", default="\t", help="field delimiter (default: tab)")
    p.add_argument("--gene-col", default="gene_id", help="gene-id column (default: gene_id)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column (default: variant_id)")
    p.add_argument("--gene-annot", default=os.environ.get("GENE_ANNOT", DEFAULT_ANNOT),
                   help="GENCODE gene table for HGNC symbols (gene_id, gene_name)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def _fmt(flux, n_snp):
    """A contributing ancestor gets a numeric flux and sign; a non-contributing one is blank."""
    if n_snp == 0:
        return "", ""
    return f"{flux:.10g}", ("+" if flux >= 0 else "-")


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    sig = load_association(args.association, d, args.gene_col, args.fdr)
    selected = load_selected(args.selected, d, args.gene_col)
    trans_w = load_trans_weights(args.weights, d, args.gene_col, args.variant_col, selected)
    per_gene, parents = load_trans_features(args.trans_features, d, args.gene_col, args.variant_col)

    wanted = {v for g in sig for _, _, v in per_gene.get(g, [])}
    sd = load_genotype_sd(args.genotype, d, args.variant_col, wanted)
    gwas_z = load_gwas_z(args.gwas, d, args.variant_col)

    label_map = gl.load_label_map(gl.default_annot(args.gene_annot), "gene_id", "gene_name")

    # One row per contributing ancestor -> disease-gene attribution, in the Supplementary Table S4
    # column layout. Each ancestor is attributed DIRECTLY to the disease gene: hop=1 if it is a
    # direct GRN parent of the gene, hop=2 otherwise (a grandparent, its hop-2 shares summed back
    # to the gene by gene_flux). An ancestor that is both a parent and a grandparent of the same
    # gene is emitted once at hop=1, so sum_A flux(A -> gene) == z_trans(gene).
    edges = []
    diseases = set()
    for gene, (sigma_g, _z_trans) in sig.items():
        rows = per_gene.get(gene)
        if not rows or sigma_g == 0.0:
            continue
        flux = gene_flux(rows, trans_w.get(gene, {}), sd, gwas_z, sigma_g)
        g_parents = parents.get(gene, set())
        for A, (fl, n_snp) in flux.items():
            if n_snp == 0:                      # no weighted, GWAS-matched SNP -> no flux delivered
                continue
            hop = 1 if A in g_parents else 2
            edges.append([gl.resolve(A, label_map), A,
                          gl.resolve(gene, label_map), gene,
                          args.tissue, str(hop), repr(float(fl))])
            diseases.add(gene)

    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"flux_{args.tissue}_{args.trait}.csv")
    header = ["regulator", "regulator_gene_id", "target", "target_gene_id", "tissue", "hop", "flux"]
    with open(out, "wt") as fh:
        fh.write(",".join(header) + "\n")
        for e in edges:
            fh.write(",".join(e) + "\n")

    n1 = sum(1 for e in edges if e[5] == "1")
    n2 = len(edges) - n1
    print(f"[build_flux_map] tissue={args.tissue} trait={args.trait} "
          f"disease_genes={len(diseases)} edges={len(edges)} (hop1={n1} hop2={n2}) -> {out}")


if __name__ == "__main__":
    main()
