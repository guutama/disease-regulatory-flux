/*
 * Optional held-out split -- per tissue, partition the samples into a train and a test fold.
 *
 * Runs only when params.train_fraction < 1 (see flux.nf). The sample ids come from the tissue's
 * hard-call genotype header (Step 5), so no separate sample list is needed. With a fraction below
 * 1 the test fold lets Step 7 score expression predictors out-of-sample; a fraction of 1 keeps the
 * whole cohort for training (this step is then skipped and Step 7 uses PSIS leave-one-out).
 *
 * Emits per tissue:
 *   <tissue>.train.samples / <tissue>.test.samples   one sample id per line
 *   <tissue>.split_qc.json                           counts + seed/fraction manifest
 *
 * Required params (see nextflow.config): train_fraction, split_seed.
 */
process SPLIT_SAMPLES {
    tag "${tissue} (train=${params.train_fraction})"
    publishDir "${params.outdir}/splits", mode: 'copy'

    cpus   1
    memory '2 GB'
    time   '20m'

    input:
    tuple val(tissue), path(genotype)

    output:
    tuple val(tissue), path("${tissue}.train.samples"), path("${tissue}.test.samples"), emit: splits
    path  "${tissue}.split_qc.json",                                                     emit: qc

    script:
    """
    build_train_test_splits.py \\
        --genotype       ${genotype} \\
        --tissue         ${tissue} \\
        --train-fraction ${params.train_fraction} \\
        --split-seed     ${params.split_seed} \\
        --outdir         .
    """
}
