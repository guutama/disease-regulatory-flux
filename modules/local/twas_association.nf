/*
 * Step 9 -- summary-statistic TWAS association with a cis/trans decomposition, per tissue.
 *
 * Takes a tissue's expression-model weights (Step 7) and hard-call genotype (Step 5) and an
 * ALT-aligned GWAS for one trait (Step 8). For each gene it forms
 * z_TWAS = sum_j w_j * sd_j * z_gwas_j / sigma_g, splits it into a cis and a trans component,
 * computes a two-sided p-value and a per-tissue BH-FDR. See bin/twas_association.py.
 *
 * Emits per (tissue, trait):
 *   association_<tissue>_<trait>.tsv   one row per gene: z_twas, z_cis, z_trans, p, p_adj, ...
 *
 * Required params (see nextflow.config): input_delimiter, gene_col.
 */
process TWAS_ASSOCIATION {
    tag "${tissue}:${trait}"
    publishDir "${params.outdir}/association", mode: 'copy'

    cpus   1
    memory '8 GB'
    time   '2h'

    input:
    tuple val(tissue), val(trait), path(weights), path(genotype), path(gwas)

    output:
    tuple val(tissue), val(trait), path("association_${tissue}_${trait}.tsv"), emit: association

    script:
    """
    twas_association.py \\
        --weights    ${weights} \\
        --genotype   ${genotype} \\
        --gwas       ${gwas} \\
        --tissue     ${tissue} \\
        --trait      ${trait} \\
        --outdir     . \\
        --delimiter  '${params.input_delimiter}' \\
        --gene-col   '${params.gene_col}'
    """
}
