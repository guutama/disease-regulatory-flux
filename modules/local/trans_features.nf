/*
 * Step 6 -- build each gene's trans-feature SNP set (the trans channel) for one tissue.
 *
 * Takes the reconstructed network (Step 3) and the LD-pruned per-gene cis-SNP set (Step 5),
 * and for each gene walks the network upward to its regulators, collecting their pruned cis
 * SNPs up to a chosen number of hops. Those SNPs are the gene's trans features. See
 * bin/trans_features.py for the contract.
 *
 * Emits per tissue:
 *   <tissue>.trans_features.tsv.gz       gene_id, hop, source_gene_id, variant_id
 *   <tissue>.trans_features_summary.tsv  per gene: number of regulators and unique SNPs
 *   <tissue>.trans_features_qc.tsv       tissue-level totals
 *
 * Required params (see nextflow.config): gene_col, trans_max_hop.
 */
process TRANS_FEATURES {
    tag "${tissue}"
    publishDir "${params.outdir}/trans_features", mode: 'copy'

    cpus   1
    memory '8 GB'
    time   '2h'

    input:
    tuple val(tissue), path(grn), path(cis_pruned)

    output:
    tuple val(tissue), path("${tissue}.trans_features.tsv.gz"),  emit: trans_features
    path  "${tissue}.trans_features_summary.tsv",                emit: summary
    path  "${tissue}.trans_features_qc.tsv",                     emit: qc

    script:
    """
    trans_features.py \\
        --grn          ${grn} \\
        --cis-pruned   ${cis_pruned} \\
        --tissue       ${tissue} \\
        --outdir       . \\
        --max-hop      ${params.trans_max_hop} \\
        --gene-col     '${params.gene_col}'
    """
}
