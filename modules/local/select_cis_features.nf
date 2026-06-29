/*
 * Step 4 -- build each gene's cis-eQTL SNP set (the cis channel) for one tissue.
 *
 * Takes the tissue's harmonised cis-eQTL table (the output of Step 1) and collects, per
 * gene, its cis-associated variants, dropping duplicate pairs. The eQTL table is already
 * cis (the eQTL analysis windowed it), so nothing is re-windowed here. If cis_max_pvalue is
 * set and the table has a p-value column, weaker associations are dropped; otherwise every
 * association is trusted. See bin/select_cis_features.py for the contract.
 *
 * Emits per tissue:
 *   <tissue>.cis_features.tsv.gz       gene_id, variant_id (one row per kept pair)
 *   <tissue>.cis_features_summary.tsv  per gene: number of cis SNPs
 *   <tissue>.cis_features_qc.tsv       pairs read, kept and dropped
 *
 * Required params (see nextflow.config): input_delimiter, gene_col; optional cis_max_pvalue.
 */
process SELECT_CIS_FEATURES {
    tag "${tissue}"
    publishDir "${params.outdir}/cis_features", mode: 'copy'

    cpus   1
    memory '4 GB'
    time   '1h'

    input:
    tuple val(tissue), path(eqtl)

    output:
    tuple val(tissue), path("${tissue}.cis_features.tsv.gz"),  emit: cis_features
    path  "${tissue}.cis_features_summary.tsv",                emit: summary
    path  "${tissue}.cis_features_qc.tsv",                     emit: qc

    script:
    def pval_opt = params.cis_max_pvalue ? "--max-pvalue ${params.cis_max_pvalue}" : ''
    """
    select_cis_features.py \\
        --eqtl       ${eqtl} \\
        --tissue     ${tissue} \\
        --outdir     . \\
        --delimiter  '${params.input_delimiter}' \\
        --gene-col   '${params.gene_col}' \\
        ${pval_opt}
    """
}
