/*
 * Step 5 -- LD-prune each gene's cis-eQTL SNP set and emit a hard-call genotype matrix for
 * one tissue.
 *
 * Takes the per-gene cis-SNP set (Step 4), the harmonised dosage matrix (Step 1) and the
 * harmonised eQTL table (Step 1, for p-value ordering). Per gene it walks the cis SNPs from
 * most to least significant and keeps a SNP only if it is not in LD (r^2 above a threshold)
 * with a SNP already kept; constant SNPs are dropped. It then writes a 0/1/2 hard-call
 * matrix for the SNPs that survive across all genes. See bin/ld_prune.py for the contract.
 *
 * Emits per tissue:
 *   <tissue>.cis_snps_pruned.tsv.gz   gene_id, variant_id (kept SNPs per gene)
 *   <tissue>.genotype_012.tsv.gz      variant_id + samples, 0/1/2 hard calls (kept SNPs)
 *   <tissue>.ld_prune_summary.tsv     per gene: cis SNPs in, kept, dropped
 *   <tissue>.ld_prune_qc.tsv          tissue-level totals
 *
 * Required params (see nextflow.config): input_delimiter, gene_col, ld_r2_threshold.
 */
process LD_PRUNE {
    tag "${tissue}"
    publishDir "${params.outdir}/ld_pruned", mode: 'copy'

    cpus   1
    memory '8 GB'
    time   '4h'

    input:
    tuple val(tissue), path(cis_features), path(genotype), path(eqtl)

    output:
    tuple val(tissue), path("${tissue}.cis_snps_pruned.tsv.gz"),
                       path("${tissue}.genotype_012.tsv.gz"),   emit: pruned
    path  "${tissue}.ld_prune_summary.tsv",                     emit: summary
    path  "${tissue}.ld_prune_qc.tsv",                          emit: qc

    script:
    """
    ld_prune.py \\
        --cis-features     ${cis_features} \\
        --genotype         ${genotype} \\
        --eqtl             ${eqtl} \\
        --tissue           ${tissue} \\
        --outdir           . \\
        --ld-r2-threshold  ${params.ld_r2_threshold} \\
        --delimiter        '${params.input_delimiter}' \\
        --gene-col         '${params.gene_col}'
    """
}
