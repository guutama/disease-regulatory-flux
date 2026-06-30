/*
 * Step 10 -- disease-regulatory flux map: attribute each disease gene's trans signal to its
 * upstream regulators, per tissue and trait.
 *
 * Takes the association (Step 9), the trans-channel weights (Step 7), the trans features
 * (Step 6), the hard-call genotype (Step 5) and the ALT-aligned GWAS (Step 8). For each
 * significant disease gene it decomposes z_trans into signed per-regulator fluxes
 * flux(A -> g) = (1/sigma_g) sum_{s in E(A) cap T(g)} w*sd*z (a SNP's term split evenly across
 * the ancestors it instruments, so the fluxes sum to z_trans(g)), and lays them out as a
 * path-preserving network: parent -> disease_gene (hop 1) and grandparent -> intermediate_parent
 * (hop 2). See bin/build_flux_map.py.
 *
 * Emits per (tissue, trait):
 *   flux_<tissue>_<trait>.tsv   one row per edge: source_gene, target_gene, disease_gene, hop,
 *                               n_snp, flux, sign, tissue, trait
 *
 * Required params (see nextflow.config): input_delimiter, gene_col, assoc_fdr.
 */
process BUILD_FLUX_MAP {
    tag "${tissue}:${trait}"
    publishDir "${params.outdir}/flux_map", mode: 'copy'

    cpus   1
    memory '8 GB'
    time   '2h'

    input:
    tuple val(tissue), val(trait), path(association), path(weights), path(trans_features),
          path(genotype), path(gwas)

    output:
    tuple val(tissue), val(trait), path("flux_${tissue}_${trait}.tsv"), emit: flux

    script:
    """
    build_flux_map.py \\
        --association     ${association} \\
        --weights         ${weights} \\
        --trans-features  ${trans_features} \\
        --genotype        ${genotype} \\
        --gwas            ${gwas} \\
        --tissue          ${tissue} \\
        --trait           ${trait} \\
        --outdir          . \\
        --fdr             ${params.assoc_fdr} \\
        --delimiter       '${params.input_delimiter}' \\
        --gene-col        '${params.gene_col}'
    """
}
