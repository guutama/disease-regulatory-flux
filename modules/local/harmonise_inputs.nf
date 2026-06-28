/*
 * Step 1 — harmonise the three reference inputs of one tissue.
 *
 * Takes the tissue's expression matrix, genotype dosage matrix and cis-eQTL map and aligns
 * them onto one shared set of samples: all expressed genes are kept (every gene is a
 * potential network target), the genotype is reduced to the variants present in the eQTL
 * map, and expression and genotype are restricted to the samples they share, in one
 * consistent order. See bin/harmonise_inputs.py for the contract and the variant-key rule
 * (rsID if present in both files, else chr:pos:ref:alt).
 *
 * Emits per tissue:
 *   <tissue>.expression.tsv.gz   all genes, shared samples
 *   <tissue>.genotype.tsv.gz     eQTL variants, shared samples
 *   <tissue>.eqtl.tsv.gz         cis pairs matched to genotype and expression
 *   <tissue>.harmonise_qc.tsv    kept/dropped counts for samples, genes and variants
 *
 * Required params (see nextflow.config): input_delimiter, gene_col, rsid_col, chrom_col,
 * pos_col, ref_col, alt_col, eqtl_beta_col, eqtl_se_col, eqtl_p_col.
 */
process HARMONISE_INPUTS {
    tag "${tissue}"
    publishDir "${params.outdir}/harmonised", mode: 'copy'

    cpus   1
    memory '8 GB'
    time   '4h'

    input:
    tuple val(tissue), path(expression), path(genotype), path(eqtl)

    output:
    tuple val(tissue), path("${tissue}.expression.tsv.gz"),
                       path("${tissue}.genotype.tsv.gz"),
                       path("${tissue}.eqtl.tsv.gz"),    emit: aligned
    path  "${tissue}.harmonise_qc.tsv",                  emit: qc

    script:
    """
    harmonise_inputs.py \\
        --expression      ${expression} \\
        --genotype        ${genotype} \\
        --eqtl            ${eqtl} \\
        --out-expression  ${tissue}.expression.tsv.gz \\
        --out-genotype    ${tissue}.genotype.tsv.gz \\
        --out-eqtl        ${tissue}.eqtl.tsv.gz \\
        --qc              ${tissue}.harmonise_qc.tsv \\
        --delimiter       '${params.input_delimiter}' \\
        --gene-col        '${params.gene_col}' \\
        --rsid-col        '${params.rsid_col}' \\
        --chrom-col       '${params.chrom_col}' \\
        --pos-col         '${params.pos_col}' \\
        --ref-col         '${params.ref_col}' \\
        --alt-col         '${params.alt_col}' \\
        --beta-col        '${params.eqtl_beta_col}' \\
        --se-col          '${params.eqtl_se_col}' \\
        --p-col           '${params.eqtl_p_col}'
    """
}
