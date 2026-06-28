#!/usr/bin/env nextflow
/*
 * Disease-regulatory flux map -- pipeline entrypoint.
 *
 * Stages:
 *   Step 1  HARMONISE_INPUTS  -- align expression, genotype and eQTL per tissue.
 *
 * Input: a samplesheet CSV (--samplesheet) with a header and one row per tissue:
 *
 *     tissue,expression,genotype,eqtl
 *     AOR,/data/AOR.expr.tsv.gz,/data/AOR.geno.tsv.gz,/data/AOR.eqtl.tsv.gz
 *
 * where the three paths are that tissue's reference files (see README).
 *
 * Run:
 *     nextflow run flux.nf --samplesheet samples.csv --outdir results_flux
 */
nextflow.enable.dsl = 2

include { HARMONISE_INPUTS } from './modules/local/harmonise_inputs.nf'

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
}
