/*
 * Step 7 -- fit a channel-aware Bayesian expression model for every gene in one tissue.
 *
 * For each gene, predicts its expression from its own cis SNPs (Step 5) and its upstream
 * regulators' cis SNPs (Step 6, the trans channel), using the harmonised expression (Step 1)
 * as the response and the 0/1/2 hard-call genotype (Step 5) as the predictors. The model is
 * fit in up to three configurations (cis_only, trans_only, cis_trans), each by MCMC (NumPyro)
 * and scored by PSIS leave-one-out (ArviZ); the best config by elpd is kept and a gene is
 * predictable when its selected leave-one-out R^2 reaches the threshold. The prior is set by
 * expr_model (bayes_ridge default, horseshoe or bslmm). See bin/fit_expression_models.py.
 *
 * Emits per tissue:
 *   <tissue>.expr_model_metrics.tsv.gz  per gene: per-config LOO-R^2, evidence, selected config
 *   <tissue>.expr_model_weights.tsv.gz  posterior-mean raw-dosage weights (long)
 *   <tissue>.expr_model_selected.tsv    the predictable genes
 *   <tissue>.expr_model_stats.tsv       tissue-level counts and mean LOO-R^2 per channel
 *
 * Required params (see nextflow.config): expr_model, expr_min_cis, expr_max_trans, gene_col.
 */
process FIT_EXPRESSION_MODELS {
    tag "${tissue}"
    publishDir "${params.outdir}/expression_models", mode: 'copy'

    cpus   1
    memory '16 GB'
    time   '8h'

    input:
    tuple val(tissue), path(expression), path(genotype), path(cis_pruned), path(trans_features)

    output:
    tuple val(tissue), path("${tissue}.expr_model_metrics.tsv.gz"),
                       path("${tissue}.expr_model_weights.tsv.gz"),  emit: models
    tuple val(tissue), path("${tissue}.expr_model_selected.tsv"),    emit: selected
    path  "${tissue}.expr_model_stats.tsv",                          emit: stats

    script:
    """
    fit_expression_models.py \\
        --expression      ${expression} \\
        --genotype        ${genotype} \\
        --cis-pruned      ${cis_pruned} \\
        --trans-features  ${trans_features} \\
        --tissue          ${tissue} \\
        --outdir          . \\
        --method          ${params.expr_model} \\
        --min-cis         ${params.expr_min_cis} \\
        --max-trans       ${params.expr_max_trans} \\
        --mcmc-warmup     ${params.expr_mcmc_warmup} \\
        --mcmc-samples    ${params.expr_mcmc_samples} \\
        --seed            ${params.expr_seed} \\
        --gene-col        '${params.gene_col}'
    """
}
