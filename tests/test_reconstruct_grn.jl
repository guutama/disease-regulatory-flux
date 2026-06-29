#=
Integration test for bin/reconstruct_grn.jl (GRN reconstruction with BioFindr).

A small simulated trio with known causal edges is generated in a temporary directory and
run through the real BioFindr pipeline. The test checks our contract -- that the expected
output files are written with the right columns, that a planted edge is recovered, that the
QC report is populated, and that mismatched sample counts are rejected -- not BioFindr's
inference, which is the library's own responsibility.

It uses only what reconstruct_grn.jl already loads (BioFindr, CSV, DataFrames) plus a tiny
self-contained random generator, so no extra packages need precompiling. Run from the
repository root:

    LD_LIBRARY_PATH=$JULIA_ROOT/lib:$JULIA_ROOT/lib/julia \
    JULIA_DEPOT_PATH=julia_project/depot:$HOME/.julia \
        julia.bin --project=julia_project tests/test_reconstruct_grn.jl
=#
include(joinpath(@__DIR__, "..", "bin", "reconstruct_grn.jl"))

# A small linear-congruential generator gives deterministic data without loading Random.
mutable struct LCG; s::UInt64; end
function unif(r::LCG)
    r.s = 6364136223846793005 * r.s + 1442695040888963407
    return (r.s >> 11) / Float64(UInt64(1) << 53)
end
gauss(r::LCG) = sqrt(-2 * log(unif(r) + 1e-12)) * cos(6.283185307179586 * unif(r))
dosage(r::LCG) = Int(floor(3 * unif(r)))            # 0, 1 or 2

"""Write a synthetic dX/dG/dE trio for tissue `T` into `dir`. Three regulators are each
driven by their instrument SNP; four targets are driven by a regulator (the planted edges
G1->G4, G1->G5, G2->G6, G3->G7), and G8 is pure noise."""
function write_trio(dir::AbstractString; n::Int=60)
    r = LCG(UInt64(7))
    genes = ["G$i" for i in 1:8]
    instr = ["G1" => "rsA", "G2" => "rsB", "G3" => "rsC"]
    geno = Dict(v => [dosage(r) for _ in 1:n] for (_, v) in instr)
    e = Dict{String,Vector{Float64}}()
    e["G1"] = [1.2 * geno["rsA"][i] + gauss(r) for i in 1:n]
    e["G2"] = [1.1 * geno["rsB"][i] + gauss(r) for i in 1:n]
    e["G3"] = [1.0 * geno["rsC"][i] + gauss(r) for i in 1:n]
    e["G4"] = [0.9 * e["G1"][i] + gauss(r) for i in 1:n]
    e["G5"] = [0.8 * e["G1"][i] + gauss(r) for i in 1:n]
    e["G6"] = [0.9 * e["G2"][i] + gauss(r) for i in 1:n]
    e["G7"] = [0.9 * e["G3"][i] + gauss(r) for i in 1:n]
    e["G8"] = [gauss(r) for _ in 1:n]
    CSV.write(joinpath(dir, "T.dX.csv"), DataFrame((g => e[g] for g in genes)...))
    CSV.write(joinpath(dir, "T.dG.csv"), DataFrame((v => geno[v] for (_, v) in instr)...))
    CSV.write(joinpath(dir, "T.dE.csv"),
              DataFrame(variant_id = [v for (_, v) in instr], gene_id = [g for (g, _) in instr]))
end

function main()
    mktempdir() do dir
        out = joinpath(dir, "out")
        write_trio(dir)
        reconstruct_grn("T", dir, out; combination="orig", fdr=0.15, method="kde")

        edges_file = joinpath(out, "T.grn_edges.csv")
        @assert isfile(edges_file) "missing $edges_file"
        @assert isfile(joinpath(out, "T.grn_edges_all.csv.gz")) "missing master edge table"
        @assert isfile(joinpath(out, "T.grn_qc.tsv")) "missing QC report"

        edges = CSV.read(edges_file, DataFrame)
        for col in ("Source", "Target", "Probability", "qvalue", "inDAG_greedy_edges")
            @assert col in names(edges) "edges file missing column $col"
        end
        @assert nrow(edges) > 0 "no edges written"
        @assert all(edges.inDAG_greedy_edges) "non-DAG edge written"
        valid = Set("G$i" for i in 1:8)
        @assert all(in(valid), edges.Source) "unexpected Source gene"
        @assert all(in(valid), edges.Target) "unexpected Target gene"

        recovered = Set((row.Source, row.Target) for row in eachrow(edges))
        @assert ("G2", "G6") in recovered "planted edge G2->G6 not recovered"

        qc = Dict(r[1] => r[2] for r in
                  CSV.File(joinpath(out, "T.grn_qc.tsv"); delim='\t', header=["metric", "value"]))
        @assert qc["combination"] == "orig" "QC combination wrong"
        @assert parse(Int, qc["regulators"]) == 3 "QC regulators != 3"
        @assert parse(Int, qc["nodes"]) >= 1 "QC nodes empty"
    end

    # mismatched sample counts between dX and dG are rejected
    threw = false
    mktempdir() do dir
        CSV.write(joinpath(dir, "T.dX.csv"), DataFrame(G1=[0.1, 0.2, 0.3], G2=[1.0, 1.1, 1.2]))
        CSV.write(joinpath(dir, "T.dG.csv"), DataFrame(rsA=[0, 1]))      # 2 rows, not 3
        CSV.write(joinpath(dir, "T.dE.csv"), DataFrame(variant_id=["rsA"], gene_id=["G1"]))
        try
            reconstruct_grn("T", dir, joinpath(dir, "out"))
        catch
            threw = true
        end
    end
    @assert threw "mismatched dX/dG sample counts were not rejected"

    println("ALL TESTS PASSED")
end

main()
