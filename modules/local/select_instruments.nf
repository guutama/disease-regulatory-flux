/*
 * Step 2 -- select each regulator's genetic instrument and build the GRN inputs for one
 * tissue.
 *
 * Takes the tissue's harmonised trio (the output of Step 1) and, for every gene with a
 * cis-eQTL, picks its lead SNP (smallest p-value; ties broken by largest effect size then
 * variant id) as the instrument. It then writes the three matrices the GRN reconstruction
 * step (findr) reads: the instrument list (dE), the instruments' genotype (dG, samples x
 * SNPs) and the expression of all genes (dX, samples x genes), with dX and dG sharing one
 * sample order. See bin/select_instruments.py for the contract.
 *
 * Emits per tissue:
 *   <tissue>.dX.csv             expression, samples x all genes
 *   <tissue>.dG.csv             instrument genotype, samples x SNPs
 *   <tissue>.dE.csv             instrument list, variant_id,gene_id pairs
 *   <tissue>.instruments_qc.tsv samples, genes and instruments selected
 *
 * Required params (see nextflow.config): input_delimiter, gene_col.
 */
process SELECT_INSTRUMENTS {
    tag "${tissue}"
    publishDir "${params.outdir}/grn_inputs", mode: 'copy'

    cpus   1
    memory '8 GB'
    time   '4h'

    input:
    tuple val(tissue), path(expression), path(genotype), path(eqtl)

    output:
    tuple val(tissue), path("${tissue}.dX.csv"),
                       path("${tissue}.dG.csv"),
                       path("${tissue}.dE.csv"),    emit: grn_inputs
    path  "${tissue}.instruments_qc.tsv",           emit: qc

    script:
    """
    select_instruments.py \\
        --expression  ${expression} \\
        --genotype    ${genotype} \\
        --eqtl        ${eqtl} \\
        --out-dx      ${tissue}.dX.csv \\
        --out-dg      ${tissue}.dG.csv \\
        --out-de      ${tissue}.dE.csv \\
        --qc          ${tissue}.instruments_qc.tsv \\
        --delimiter   '${params.input_delimiter}' \\
        --gene-col    '${params.gene_col}'
    """
}
