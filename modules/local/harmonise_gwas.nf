/*
 * Step 8 -- harmonise a GWAS summary-statistics file to the pipeline's ALT allele.
 *
 * Takes a raw GWAS file for one trait and the pipeline's variant universe (every tissue's
 * genotype_012 matrix from Step 5, which carries chromosome, position, ref and alt). It
 * matches each variant by chromosome and position, aligns the GWAS effect (and z) to the
 * ALT allele, drops strand-ambiguous and unmatched variants, and writes one ALT-aligned
 * record per variant. See bin/harmonise_gwas.py for the contract.
 *
 * Emits per trait:
 *   <trait>.gwas.tsv.gz    variant_id, z, beta, se, p_value, n  (z and beta on the ALT scale)
 *   <trait>.gwas_qc.tsv    trait-level totals (aligned / dropped)
 *
 * GWAS column names default to the GWAS-Catalog harmonised layout (chromosome, position,
 * effect_allele, other_allele, beta, se, pvalue, n); override per cohort if they differ.
 *
 * Required params (see nextflow.config): input_delimiter.
 */
process HARMONISE_GWAS {
    tag "${trait}"
    publishDir "${params.outdir}/gwas", mode: 'copy'

    cpus   1
    memory '8 GB'
    time   '2h'

    input:
    tuple val(trait), path(gwas)
    path  variants

    output:
    tuple val(trait), path("${trait}.gwas.tsv.gz"), emit: harmonised
    path  "${trait}.gwas_qc.tsv",                    emit: qc

    script:
    """
    harmonise_gwas.py \\
        --gwas       ${gwas} \\
        --variants   ${variants} \\
        --trait      ${trait} \\
        --out        ${trait}.gwas.tsv.gz \\
        --qc         ${trait}.gwas_qc.tsv \\
        --delimiter  '${params.input_delimiter}'
    """
}
