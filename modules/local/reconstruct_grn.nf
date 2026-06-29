/*
 * Step 3 -- reconstruct one tissue's gene regulatory network with BioFindr.
 *
 * Takes the tissue's dX/dG/dE matrices (the output of Step 2) and runs the canonical
 * three-function BioFindr pipeline (findr -> globalfdr! -> dagfindr!) to infer a directed
 * acyclic regulatory network: every gene is a candidate target, each regulator is oriented
 * by its instrument. See bin/reconstruct_grn.jl for the contract.
 *
 * BioFindr is a Julia package, pinned in params.julia_project_dir (Julia 1.11.3). The
 * posterior loops use Julia threads, so the task runs with -t equal to its cpu count.
 *
 * Emits per tissue:
 *   <tissue>.grn_edges.csv          the network (Source, Target, Probability, qvalue, ...)
 *   <tissue>.grn_edges_all.csv.gz   every inferred edge with its q-value, before the FDR cut
 *   <tissue>.grn_qc.tsv             node, edge, regulator, target and degree counts
 *
 * Required params (see nextflow.config): julia_bin, julia_project_dir, grn_cpus,
 * grn_combination, grn_fdr, grn_findr_method.
 */
process RECONSTRUCT_GRN {
    tag "${tissue}"
    publishDir "${params.outdir}/grn", mode: 'copy'

    cpus   { params.grn_cpus as int }
    memory '64 GB'
    time   '8h'

    input:
    tuple val(tissue), path(dX), path(dG), path(dE)

    output:
    tuple val(tissue), path("${tissue}.grn_edges.csv"),
                       path("${tissue}.grn_edges_all.csv.gz"),  emit: grn
    path  "${tissue}.grn_qc.tsv",                               emit: qc

    script:
    """
    # BioFindr's Julia (julia.bin) needs its own shared libraries on LD_LIBRARY_PATH, and
    # the project's precompiled package cache on JULIA_DEPOT_PATH.
    julia_lib="\$(cd "\$(dirname '${params.julia_bin}')/../lib" && pwd)"
    export LD_LIBRARY_PATH="\${julia_lib}:\${julia_lib}/julia:\${LD_LIBRARY_PATH:-}"
    export JULIA_DEPOT_PATH='${params.julia_project_dir}'/depot:\$HOME/.julia

    '${params.julia_bin}' --project='${params.julia_project_dir}' -t ${task.cpus} \\
        ${projectDir}/bin/reconstruct_grn.jl \\
        '${tissue}' . . \\
        '${params.grn_combination}' \\
        '${params.grn_fdr}' \\
        '${params.grn_findr_method}'
    """
}
