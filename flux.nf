#!/usr/bin/env nextflow
/*
 * Disease-regulatory flux map -- pipeline entrypoint.
 *
 * Stages:
 *   Step 1  HARMONISE_INPUTS    -- align expression, genotype and eQTL per tissue.
 *   Step 2  SELECT_INSTRUMENTS  -- pick each regulator's lead cis-eQTL and build the
 *                                 GRN reconstruction inputs (dX, dG, dE) per tissue.
 *   Step 3  RECONSTRUCT_GRN     -- infer the directed gene regulatory network per tissue
 *                                 with BioFindr.
 *   Step 4  SELECT_CIS_FEATURES -- build each gene's cis-eQTL SNP set (the cis channel)
 *                                 per tissue.
 *   Step 5  LD_PRUNE            -- LD-prune each gene's cis SNPs and emit a 0/1/2 hard-call
 *                                 genotype matrix per tissue.
 *   Step 6  TRANS_FEATURES      -- collect each gene's upstream regulators' pruned cis SNPs
 *                                 (the trans channel) by walking the network per tissue.
 *   Step 7  FIT_EXPRESSION_MODELS -- fit channel-aware Bayesian expression models (cis /
 *                                 trans / cis+trans) per gene, per tissue.
 *   Step 8  HARMONISE_GWAS      -- align a GWAS to the ALT allele (optional; per trait).
 *   Step 9  TWAS_ASSOCIATION    -- summary-statistic TWAS with a cis/trans split (optional;
 *                                 per tissue and trait).
 *
 * Input: a samplesheet CSV (--samplesheet) with a header and one row per tissue:
 *
 *     tissue,expression,genotype,eqtl
 *     AOR,/data/AOR.expr.tsv.gz,/data/AOR.geno.tsv.gz,/data/AOR.eqtl.tsv.gz
 *
 * where the three paths are that tissue's reference files (see README).
 *
 * Run:
 *     nextflow run flux.nf --samplesheet samplesheet.csv --outdir results_flux
 */
nextflow.enable.dsl = 2

include { HARMONISE_INPUTS   } from './modules/local/harmonise_inputs.nf'
include { SELECT_INSTRUMENTS } from './modules/local/select_instruments.nf'
include { RECONSTRUCT_GRN    } from './modules/local/reconstruct_grn.nf'
include { SELECT_CIS_FEATURES } from './modules/local/select_cis_features.nf'
include { LD_PRUNE           } from './modules/local/ld_prune.nf'
include { TRANS_FEATURES     } from './modules/local/trans_features.nf'
include { FIT_EXPRESSION_MODELS } from './modules/local/fit_expression_models.nf'
include { HARMONISE_GWAS     } from './modules/local/harmonise_gwas.nf'
include { TWAS_ASSOCIATION   } from './modules/local/twas_association.nf'

workflow {
    if( !params.samplesheet )
        error "Provide --samplesheet : a CSV with columns tissue,expression,genotype,eqtl"

    inputs = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.tissue,
            file(row.expression, checkIfExists: true),
            file(row.genotype,   checkIfExists: true),
            file(row.eqtl,       checkIfExists: true)) }

    HARMONISE_INPUTS(inputs)
    SELECT_INSTRUMENTS(HARMONISE_INPUTS.out.aligned)
    RECONSTRUCT_GRN(SELECT_INSTRUMENTS.out.grn_inputs)

    // cis channel: build each gene's cis-eQTL SNP set from the harmonised eQTL.
    SELECT_CIS_FEATURES(HARMONISE_INPUTS.out.aligned.map { t, x, g, e -> tuple(t, e) })

    // LD-prune the cis SNP set, using the harmonised genotype and eQTL (joined by tissue).
    ld_in = SELECT_CIS_FEATURES.out.cis_features
        .join(HARMONISE_INPUTS.out.aligned)
        .map { t, cis, x, g, e -> tuple(t, cis, g, e) }
    LD_PRUNE(ld_in)

    // trans channel: walk the GRN upward to collect each gene's ancestors' pruned cis SNPs.
    tf_in = RECONSTRUCT_GRN.out.grn
        .join(LD_PRUNE.out.pruned)
        .map { t, edges, edges_all, cis_pruned, geno012 -> tuple(t, edges, cis_pruned) }
    TRANS_FEATURES(tf_in)

    // expression models: predict each gene from its cis SNPs and its regulators' cis SNPs.
    fit_in = HARMONISE_INPUTS.out.aligned.map { t, x, g, e -> tuple(t, x) }
        .join(LD_PRUNE.out.pruned)                       // t, expr, cis_pruned, geno012
        .join(TRANS_FEATURES.out.trans_features)         // t, expr, cis_pruned, geno012, trans
        .map { t, expr, cis_pruned, geno012, trans -> tuple(t, expr, geno012, cis_pruned, trans) }
    FIT_EXPRESSION_MODELS(fit_in)

    // association against GWAS traits (optional; runs only when --gwas-samplesheet is given).
    if( params.gwas_samplesheet ) {
        gwas_inputs = Channel
            .fromPath(params.gwas_samplesheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> tuple(row.trait, file(row.gwas, checkIfExists: true)) }

        // variant universe for GWAS alignment: every tissue's hard-call genotype.
        variants = LD_PRUNE.out.pruned.map { t, cis_pruned, geno012 -> geno012 }.collect()
        HARMONISE_GWAS(gwas_inputs, variants)

        // each tissue's weights + genotype, crossed with each trait's harmonised GWAS.
        weights_geno = FIT_EXPRESSION_MODELS.out.models
            .map { t, metrics, weights -> tuple(t, weights) }
            .join(LD_PRUNE.out.pruned.map { t, cis_pruned, geno012 -> tuple(t, geno012) })
        assoc_in = weights_geno.combine(HARMONISE_GWAS.out.harmonised)
            .map { t, weights, geno012, trait, gwas -> tuple(t, trait, weights, geno012, gwas) }
        TWAS_ASSOCIATION(assoc_in)
    }
}
