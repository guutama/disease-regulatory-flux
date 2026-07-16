/*
 * Step 3 -- reconstruct one tissue's gene regulatory network with findr (python libfindr).
 *
 * Takes the tissue's dX/dG/dE matrices (the output of Step 2) and runs the canonical
 * three-step findr pipeline (posteriors -> globalfdr -> dagfindr) to infer a directed
 * acyclic regulatory network: every gene is a candidate target, each regulator is oriented
 * by its instrument. See bin/run_findr_py.py for the contract.
 *
 * findr calls Lingfei Wang's libfindr C core (params.findr_libpath), which links against
 * GSL (params.findr_gsl_dir, exported so the shared library loads). The posterior step is
 * threaded with --nth equal to the task's cpu count.
 *
 * Emits per tissue:
 *   <tissue>.grn_edges.csv          the network (Source, Target, Probability, qvalue, ...)
 *   <tissue>.grn_edges_all.csv.gz   every inferred edge with its q-value, before the FDR cut
 *   <tissue>.grn_qc.tsv             node, edge, regulator, target and degree counts
 *
 * Required params (see nextflow.config): findr_libpath, findr_gsl_dir, grn_cpus,
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
    # libfindr.so links against GSL; expose its directory so the shared library loads.
    export FINDR_LIBPATH='${params.findr_libpath}'
    export GSL_LIB_DIR='${params.findr_gsl_dir}'
    export LD_LIBRARY_PATH="\${GSL_LIB_DIR}:\${LD_LIBRARY_PATH:-}"

    run_findr_py.py \\
        '${tissue}' . . \\
        '${params.grn_combination}' \\
        '${params.grn_fdr}' \\
        '${params.grn_findr_method}' \\
        --nth ${task.cpus} \\
        --libpath '${params.findr_libpath}'
    """
}
