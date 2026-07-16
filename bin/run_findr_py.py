#!/usr/bin/env python3.9
"""
Reconstruct one tissue's gene regulatory network with findr (the python libfindr package).

Reads the three matrices written by the instrument-selection step (comma-separated):

    dX  expression, samples x genes         (column names = gene ids)
    dG  instrument genotype, samples x SNPs  (column names = variant ids)
    dE  instrument list, variant_id,gene_id  (the SNP in column 1, the gene in column 2)

and runs the canonical three-step findr pipeline:

    posteriors   for every (regulator, target) pair, the posterior probability of a directed
                 regulatory edge, using the regulator's instrument to orient the edge;
    globalfdr    turn the posteriors into q-values across all non-self edges and keep the
                 edges at the requested false-discovery rate;
    dagfindr     remove cycles greedily (add edges by decreasing probability, skip any that
                 would close a cycle) to leave a directed acyclic network.

Every gene in dX is a candidate target; the regulators are the genes listed in dE, so each
edge starts from an instrument-anchored source. The libfindr C core supernormalises the
expression internally (rank-based inverse-normal), so raw expression is passed unchanged.

The posterior step is the expensive one and runs once. Its full edge table (every edge with
its q-value, before the FDR cut and cycle removal) is written as well, so a network at a
different FDR can be derived later without recomputing the posteriors.

Outputs (in <out_dir>):
    <tissue>.grn_edges.csv         the network: the directed acyclic edges kept at the
                                   requested FDR (Source, Target, Probability, qvalue,
                                   Source_idx, Target_idx, inDAG_greedy_edges);
    <tissue>.grn_edges_all.csv.gz  every inferred edge with its q-value, before the FDR cut,
                                   in decreasing probability;
    <tissue>.grn_qc.tsv            counts of nodes, edges, regulators and targets, and the
                                   mean and maximum out- and in-degree.

Usage:
    run_findr_py.py <tissue> <input_dir> <out_dir> \
        [<combination=orig>] [<fdr=0.15>] [<method=kde>] \
        [--nth N] [--na K] [--libpath PATH] [--no-all-edges]

  <combination> is the findr test combination that builds the edge posterior from the five
  sub-tests: "orig" = 0.5*(P2*P5 + P4), "IV" = P2*P5, "mediation" = P2*P3.
  <method> is a provenance label recorded in the QC report; libfindr's posterior conversion
  is fixed and the label does not change it.

The libfindr shared library is located via --libpath or the FINDR_LIBPATH environment
variable. Where GSL is not already on LD_LIBRARY_PATH, its library directory can be given in
the GSL_LIB_DIR environment variable so libfindr.so still links.
"""
import argparse
import ctypes
import glob
import gzip
import os
import sys
import time

import numpy as np
import pandas as pd

DEFAULT_LIBPATH = os.environ.get("FINDR_LIBPATH", "/cluster/projects/nn1015k/findr/libfindr.so")


def gsl_lib_dirs():
    """Directories to search for the GSL shared libraries: GSL_LIB_DIR first, then the
    common cluster module locations (globs allowed)."""
    dirs = []
    env = os.environ.get("GSL_LIB_DIR")
    if env:
        dirs.append(env)
    dirs += ["/cluster/software/GSL/2.7-GCC-13.2.0/lib",
             "/cluster/software/GSL/2.7-*/lib",
             "/cluster/software/GSL/*/lib"]
    return dirs


def preload_gsl():
    """Best-effort preload of libgslcblas / libgsl so libfindr.so links even when the GSL
    module is not loaded (LD_LIBRARY_PATH unset)."""
    for name in ("libgslcblas.so.0", "libgsl.so.25"):
        try:
            ctypes.CDLL(name, mode=ctypes.RTLD_GLOBAL)
            continue
        except OSError:
            pass
        loaded = False
        for pat in gsl_lib_dirs():
            hits = sorted(glob.glob(os.path.join(pat, name)))
            if hits:
                ctypes.CDLL(hits[0], mode=ctypes.RTLD_GLOBAL)
                loaded = True
                break
        if not loaded:
            sys.stderr.write(
                "WARNING: could not preload {}; relying on LD_LIBRARY_PATH "
                "(set GSL_LIB_DIR or module load GSL).\n".format(name))


def init_findr(libpath, nth):
    import findr
    return findr.lib(path=libpath, loglv=6, rs=0, nth=nth)


def getpairs(genes_X, snps_G, dE):
    """(colG, colX) 0-based index pairs, one per dE row whose gene is in dX, sorted by the
    source-gene column colX. colG indexes the instrument SNP in dG, colX its regulator gene
    in dX."""
    gene2idx = {g: i for i, g in enumerate(genes_X)}
    snp2idx = {s: i for i, s in enumerate(snps_G)}
    idG = dE.iloc[:, 0].tolist()
    idX = dE.iloc[:, 1].tolist()
    pairs = []
    n_missing_gene = 0
    for sG, sX in zip(idG, idX):
        jX = gene2idx.get(sX)
        if jX is None:
            n_missing_gene += 1
            continue
        jG = snp2idx.get(sG)
        if jG is None:
            raise ValueError('Variant ID "{}" not found in dG columns.'.format(sG))
        pairs.append((jG, jX))
    if n_missing_gene:
        sys.stderr.write("  {} eQTL rows dropped (gene not in dX)\n".format(n_missing_gene))
    pairs = np.array(pairs, dtype=np.int64)
    pairs = pairs[np.argsort(pairs[:, 1], kind="stable")]
    return pairs


def qvalue(P):
    """Storey-style q-values for posterior probabilities P: sort descending, take the
    cumulative mean of (1 - P), enforce a monotone non-decreasing tail, clamp to [0, 1] and
    restore the original order."""
    P = np.asarray(P, dtype=np.float64)
    n = P.size
    order = np.argsort(-P, kind="stable")
    Psorted = P[order]
    csum = np.cumsum(Psorted)
    qsort = 1.0 - csum / np.arange(1, n + 1)
    qsort = np.minimum.accumulate(qsort[::-1])[::-1]
    np.clip(qsort, 0.0, 1.0, out=qsort)
    qval = np.empty(n, dtype=np.float64)
    qval[order] = qsort
    return qval


def dagfindr_greedy_edges(df):
    """Greedy-edges DAG: add edges by decreasing probability and drop any that would close a
    cycle (only edges whose target can itself be a source can create one). Mutates+returns df
    with Source_idx / Target_idx / inDAG_greedy_edges columns, plus the vertex map name2idx
    (1-based) and the resulting acyclic graph G."""
    import networkx as nx
    df = df.sort_values(["qvalue", "Probability"], ascending=[True, False],
                        kind="stable").reset_index(drop=True)
    source_names = set(df["Source"].unique())
    vnames = pd.unique(pd.concat([df["Source"], df["Target"]], ignore_index=True))
    name2idx = {name: i + 1 for i, name in enumerate(vnames)}
    df["Source_idx"] = df["Source"].map(name2idx).astype(np.int64)
    df["Target_idx"] = df["Target"].map(name2idx).astype(np.int64)

    G = nx.DiGraph()
    G.add_nodes_from(name2idx.values())
    indag = np.ones(len(df), dtype=bool)
    src = df["Source_idx"].to_numpy()
    tgt = df["Target_idx"].to_numpy()
    tgt_name = df["Target"].to_numpy()
    for k in range(len(df)):
        u = int(src[k]); v = int(tgt[k])
        if tgt_name[k] in source_names:
            if u == v or nx.has_path(G, v, u):
                indag[k] = False
            else:
                G.add_edge(u, v)
        else:
            G.add_edge(u, v)
    df["inDAG_greedy_edges"] = indag
    return df, name2idx, G


def combine(res, combination):
    """Build the edge posterior from libfindr's five sub-test probabilities."""
    p2 = res["p2"]; p3 = res["p3"]; p4 = res["p4"]; p5 = res["p5"]
    if combination == "orig":
        return 0.5 * (p2 * p5 + p4)
    elif combination == "IV":
        return p2 * p5
    elif combination == "mediation":
        return p2 * p3
    raise ValueError("combination must be orig, IV, or mediation")


def write_all_edges(path_gz, PP, qfull, valid, genes, colX, chunk=5_000_000):
    """Write every non-self edge (Source, Target, Probability, qvalue) gzip-compressed, in
    decreasing probability. Rows index sources (their gene is genes[colX[row]]); columns
    index target genes. Written in chunks so the full edge list is never materialised at
    once."""
    genes_arr = np.asarray(genes, dtype=object)
    vr, vc = np.where(valid)
    prob = PP[vr, vc]
    qv = qfull[vr, vc]
    order = np.argsort(-prob, kind="stable")
    with gzip.open(path_gz, "wt") as fh:
        fh.write("Source,Target,Probability,qvalue\n")
        for s in range(0, order.size, chunk):
            idx = order[s:s + chunk]
            r = vr[idx]; c = vc[idx]
            pd.DataFrame({
                "Source": genes_arr[colX[r]],
                "Target": genes_arr[c],
                "Probability": prob[idx],
                "qvalue": qv[idx],
            }).to_csv(fh, index=False, header=False)
    return int(order.size)


def reconstruct_grn(tissue, input_dir, out_dir, combination, fdr, method,
                    nth, na, libpath, write_all=True):
    t0 = time.time()
    os.makedirs(out_dir, exist_ok=True)
    preload_gsl()
    lib = init_findr(libpath, nth)

    dX = pd.read_csv(os.path.join(input_dir, "{}.dX.csv".format(tissue)))
    dG = pd.read_csv(os.path.join(input_dir, "{}.dG.csv".format(tissue)))
    dE = pd.read_csv(os.path.join(input_dir, "{}.dE.csv".format(tissue)))
    genes = list(dX.columns)
    snps = list(dG.columns)
    if dX.shape[0] != dG.shape[0]:
        raise ValueError("dX has {} sample rows but dG has {}; their sample order must "
                         "agree.".format(dX.shape[0], dG.shape[0]))

    pairs = getpairs(genes, snps, dE)
    npairs = pairs.shape[0]
    print("[{}] dX={}x{}  dG={}x{}  regulators={}  dE={}  comb={}  fdr={}  method={}".format(
        tissue, dX.shape[0], dX.shape[1], dG.shape[0], dG.shape[1], npairs, dE.shape[0],
        combination, fdr, method))

    X = dX.to_numpy(dtype=np.float32).T           # (ngenes, ns)
    Gmat = dG.to_numpy().T                         # (nsnps, ns)
    ns = X.shape[1]; ngenes = X.shape[0]
    colG = pairs[:, 0]; colX = pairs[:, 1]

    # Pass every gene as a target, but order dt2 so the sources occupy its first `npairs`
    # rows aligned to dt, and set nodiag=True so libfindr skips each source's own row (the
    # source-vs-self diagonal). Feeding a source its own identical row instead yields a
    # degenerate LLR that corrupts that source's posterior fit. Sources are unique (one lead
    # eQTL per regulator), so the alignment is a clean bijection.
    src_in_target = set(colX.tolist())
    others = np.array([j for j in range(ngenes) if j not in src_in_target], dtype=np.int64)
    dt2_order = np.concatenate([colX, others])     # first npairs rows = sources
    dg = np.ascontiguousarray(Gmat[colG, :]).astype("u1")   # (npairs, ns) instrument dosage
    dt = np.ascontiguousarray(X[colX, :])                   # (npairs, ns) source expression
    dt2 = np.ascontiguousarray(X[dt2_order, :])             # (ngenes, ns) target expression

    print("[{}] pijs_gassist: {} sources x {} targets x {} samples ...".format(
        tissue, npairs, ngenes, ns))
    tt = time.time()
    res = lib.pijs_gassist(dg, dt, dt2, na=na, nodiag=True, memlimit=-1)
    print("  pijs_gassist ret={} in {:.1f}s".format(res["ret"], time.time() - tt))
    if res["ret"] != 0:
        sys.stderr.write("  WARNING: pijs_gassist ret={} (nonzero); validating outputs are "
                         "finite.\n".format(res["ret"]))

    # combine the sub-tests, then remap columns from the dt2_order layout back to gene index
    PP_ord = combine(res, combination).astype(np.float64)   # (npairs, ngenes) in dt2_order
    PP = np.empty_like(PP_ord)
    PP[:, dt2_order] = PP_ord
    n_nan = int(np.isnan(PP).sum()); n_neg = int((PP < 0).sum())
    if n_nan or n_neg:
        print("  sanitising PP: {} NaN, {} negative -> 0".format(n_nan, n_neg))
    PP = np.nan_to_num(PP, nan=0.0, posinf=1.0, neginf=0.0)
    np.clip(PP, 0.0, 1.0, out=PP)

    # q-values over all non-self edges (one self entry per row, at colX[k])
    self_mask = np.zeros(PP.shape, dtype=bool)
    self_mask[np.arange(npairs), colX] = True
    valid = ~self_mask
    qflat = qvalue(PP[valid])
    qfull = np.empty(PP.shape, dtype=np.float64)
    qfull[valid] = qflat
    qfull[self_mask] = np.inf

    # (a) master edge table: every non-self edge with its q-value, before the FDR cut
    all_path = os.path.join(out_dir, "{}.grn_edges_all.csv.gz".format(tissue))
    if write_all:
        n_all = write_all_edges(all_path, PP, qfull, valid, genes, colX)
        print("[{}] wrote {}  ({} edges)".format(tissue, all_path, n_all))

    # (b) keep the edges passing the FDR threshold, then greedy-edges DAG
    keep = qfull <= fdr
    rows, cols = np.where(keep)
    genes_arr = np.array(genes, dtype=object)
    dP = pd.DataFrame({
        "Source": genes_arr[colX[rows]],
        "Target": genes_arr[cols],
        "Probability": PP[rows, cols],
        "qvalue": qfull[rows, cols],
    })
    dP, name2idx, G = dagfindr_greedy_edges(dP)
    dag = dP[dP["inDAG_greedy_edges"]].copy()
    edges_path = os.path.join(out_dir, "{}.grn_edges.csv".format(tissue))
    dag[["Source", "Target", "Probability", "qvalue",
         "Source_idx", "Target_idx", "inDAG_greedy_edges"]].to_csv(edges_path, index=False)
    print("[{}] wrote {}  ({} DAG edges, {} dropped as cycles)".format(
        tissue, edges_path, len(dag), len(dP) - len(dag)))

    # (c) QC: degree summary of the reconstructed acyclic network
    outdeg = dict(G.out_degree()); indeg = dict(G.in_degree())
    nreg = sum(1 for d in outdeg.values() if d > 0)
    ntar = sum(1 for d in indeg.values() if d > 0)
    ne_g = G.number_of_edges()
    qc_path = os.path.join(out_dir, "{}.grn_qc.tsv".format(tissue))
    with open(qc_path, "w") as io:
        io.write("metric\tvalue\n")
        io.write("combination\t{}\n".format(combination))
        io.write("fdr\t{}\n".format(fdr))
        io.write("method\t{}\n".format(method))
        io.write("nodes\t{}\n".format(G.number_of_nodes()))
        io.write("edges\t{}\n".format(ne_g))
        io.write("regulators\t{}\n".format(nreg))
        io.write("targets\t{}\n".format(ntar))
        io.write("mean_out_degree\t{}\n".format(ne_g / nreg if nreg else 0.0))
        io.write("max_out_degree\t{}\n".format(max(outdeg.values(), default=0)))
        io.write("mean_in_degree\t{}\n".format(ne_g / ntar if ntar else 0.0))
        io.write("max_in_degree\t{}\n".format(max(indeg.values(), default=0)))
    print("[{}] nodes={} edges={} regulators={} targets={} -> {}  ({:.1f}s total)".format(
        tissue, G.number_of_nodes(), ne_g, nreg, ntar, qc_path, time.time() - t0))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tissue")
    ap.add_argument("input_dir", help="directory with <tissue>.dX.csv/.dG.csv/.dE.csv")
    ap.add_argument("out_dir")
    ap.add_argument("combination", nargs="?", default="orig", help="orig | IV | mediation")
    ap.add_argument("fdr", nargs="?", type=float, default=0.15,
                    help="false-discovery-rate threshold (default: 0.15)")
    ap.add_argument("method", nargs="?", default="kde",
                    help="provenance label recorded in QC (libfindr conversion is fixed)")
    ap.add_argument("--nth", type=int, default=4, help="libfindr worker threads")
    ap.add_argument("--na", type=int, default=None,
                    help="number of alleles (nvx=na+1); default: auto (max genotype)")
    ap.add_argument("--libpath", default=DEFAULT_LIBPATH,
                    help="path to libfindr.so (or set FINDR_LIBPATH)")
    ap.add_argument("--no-all-edges", action="store_true",
                    help="skip the full <tissue>.grn_edges_all.csv.gz master table")
    args = ap.parse_args()

    reconstruct_grn(args.tissue, args.input_dir, args.out_dir, args.combination,
                    args.fdr, args.method, args.nth, args.na, args.libpath,
                    write_all=not args.no_all_edges)


if __name__ == "__main__":
    main()
