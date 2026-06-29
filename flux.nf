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
}
