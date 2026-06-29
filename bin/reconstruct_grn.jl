#!/usr/bin/env julia
"""
Reconstruct one tissue's gene regulatory network with BioFindr.

Reads the three matrices written by the instrument-selection step (comma-separated):

    dX  expression, samples x genes        (column names = gene ids)
    dG  instrument genotype, samples x SNPs (column names = variant ids)
    dE  instrument list, variant_id,gene_id (the SNP in column 1, the gene in column 2)

and runs the canonical three-function BioFindr pipeline:

    dP          = findr(dX, dG, dE; colX=2, colG=1, method, combination, FDR=1.0, sorted=true)
    globalfdr!(dP; FDR, sorted=true)
    G, name2idx = dagfindr!(dP; method="greedy edges")

findr infers, for every (regulator, target) pair, the posterior probability of a directed
regulatory edge, using the regulator's instrument to orient the edge. Every gene in dX is a
candidate target; the regulators are the genes listed in dE, so each edge starts from an
instrument-anchored source. globalfdr! turns the posteriors into q-values and keeps the
edges at the requested false-discovery rate, and dagfindr! removes cycles greedily to leave
a directed acyclic network.

The posterior step is the expensive one and runs once. Its full edge table (every edge with
its q-value, before the FDR cut and cycle removal) is written as well, so a network at a
different FDR can be derived later without recomputing the posteriors.

Outputs (in <out_dir>):
    <tissue>.grn_edges.csv       the network: Source, Target and the per-edge columns from
                                 BioFindr, restricted to the directed acyclic edges kept at
                                 the requested FDR;
    <tissue>.grn_edges_all.csv.gz  every inferred edge with its q-value, before the FDR cut;
    <tissue>.grn_qc.tsv          counts of nodes, edges, regulators and targets, and the
                                 mean and maximum out- and in-degree.

Usage:
    julia --project=<julia_project> [-t N] reconstruct_grn.jl \\
        <tissue> <inputs_dir> <out_dir> [<combination=orig>] [<fdr=0.15>] [<method=kde>]

  <combination> is BioFindr's test combination: "orig", "IV" or "mediation".
  <method> is the LLR mixture fit: "kde" (robust) or "moments" (faster).
"""

using BioFindr, CSV, DataFrames, Graphs, Printf, Logging

"""
    reconstruct_grn(tissue, input_dir, out_dir; combination="orig", fdr=0.15, method="kde")

Reconstruct one tissue's GRN from the dX/dG/dE inputs in `input_dir` and write the network,
the master edge table and the QC report to `out_dir`. Returns the BioFindr graph `G`.
"""
function reconstruct_grn(tissue::AbstractString, input_dir::AbstractString,
                         out_dir::AbstractString;
                         combination::AbstractString="orig",
                         fdr::Real=0.15, method::AbstractString="kde")
    mkpath(out_dir)

    dX = CSV.read(joinpath(input_dir, "$(tissue).dX.csv"), DataFrame)   # samples x all genes
    dG = CSV.read(joinpath(input_dir, "$(tissue).dG.csv"), DataFrame)   # samples x instrument SNPs
    dE = CSV.read(joinpath(input_dir, "$(tissue).dE.csv"), DataFrame)   # variant_id, gene_id

    nrow(dX) == nrow(dG) || error(
        "dX has $(nrow(dX)) sample rows but dG has $(nrow(dG)); their sample order must agree.")
    regulators = unique(string.(dE[:, 2]))
    all(in(Set(names(dX))), regulators) || error(
        "some regulator genes in dE are absent from the expression matrix dX.")
    @printf("[%s] threads=%d  dX=%dx%d (all targets)  regulators=%d  dE=%d  comb=%s  fdr=%.2f  method=%s\n",
            tissue, Threads.nthreads(), nrow(dX), ncol(dX), length(regulators), nrow(dE),
            combination, fdr, method)

    # (1) findr -- pairwise causal posteriors, full table with q-values (FDR=1.0 keeps all)
    t0 = time()
    dP = with_logger(ConsoleLogger(stderr, Logging.Error)) do
        findr(dX, dG, dE; colX=2, colG=1, method=method, combination=combination,
              FDR=1.0, sorted=true)
    end
    @printf("[%s] findr(comb=%s) -> %d edges in %.0fs\n", tissue, combination, nrow(dP), time() - t0)

    # master edge table (every edge with its q-value, before the FDR cut) for later reuse
    all_file = joinpath(out_dir, "$(tissue).grn_edges_all.csv")
    CSV.write(all_file, dP)
    try
        run(`gzip -f $all_file`)
        println("[$tissue] wrote $(all_file).gz  ($(nrow(dP)) edges)")
    catch e
        println("[$tissue] WARNING: gzip failed ($e); left uncompressed $all_file")
    end

    # (2) globalfdr! -- keep edges at the requested FDR (in place)
    BioFindr.globalfdr!(dP; FDR=fdr, sorted=true)

    # (3) dagfindr! -- greedy cycle removal -> directed acyclic graph G
    G, _ = BioFindr.dagfindr!(dP; method="greedy edges")
    dag = filter(:inDAG_greedy_edges => identity, dP)

    edges_file = joinpath(out_dir, "$(tissue).grn_edges.csv")
    CSV.write(edges_file, dag)
    println("[$tissue] wrote $edges_file  ($(nrow(dag)) DAG edges)")

    # QC: degree summary of the reconstructed network
    outdeg = outdegree(G); indeg = indegree(G)
    nreg = count(>(0), outdeg); ntar = count(>(0), indeg); ne_g = ne(G)
    qc_file = joinpath(out_dir, "$(tissue).grn_qc.tsv")
    open(qc_file, "w") do io
        println(io, "metric\tvalue")
        println(io, "combination\t$combination")
        println(io, "fdr\t$fdr")
        println(io, "method\t$method")
        println(io, "nodes\t$(nv(G))")
        println(io, "edges\t$ne_g")
        println(io, "regulators\t$nreg")
        println(io, "targets\t$ntar")
        println(io, "mean_out_degree\t$(nreg > 0 ? ne_g / nreg : 0.0)")
        println(io, "max_out_degree\t$(maximum(outdeg; init=0))")
        println(io, "mean_in_degree\t$(ntar > 0 ? ne_g / ntar : 0.0)")
        println(io, "max_in_degree\t$(maximum(indeg; init=0))")
    end
    @printf("[%s] nodes=%d edges=%d regulators=%d targets=%d -> %s\n",
            tissue, nv(G), ne_g, nreg, ntar, qc_file)
    return G
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 3 || error(
        "Usage: julia reconstruct_grn.jl <tissue> <inputs_dir> <out_dir> " *
        "[<combination=orig>] [<fdr=0.15>] [<method=kde>]")
    reconstruct_grn(ARGS[1], ARGS[2], ARGS[3];
                    combination = length(ARGS) >= 4 ? ARGS[4] : "orig",
                    fdr         = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.15,
                    method      = length(ARGS) >= 6 ? ARGS[6] : "kde")
end
