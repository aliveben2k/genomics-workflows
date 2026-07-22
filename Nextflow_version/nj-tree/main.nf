#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { NJ_TREE } from './workflows/nj-tree'


workflow {
    if (!params.input) {
        error "Missing required parameter: --input"
    }
    if (!params.outdir) {
        error "Missing required parameter: --outdir"
    }
    if (params.window_size <= 0) {
        error "--window_size must be greater than zero"
    }

    ch_vcf = Channel.fromPath(params.input, checkIfExists: true)
    ch_pop = params.pop
        ? Channel.fromPath(params.pop, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/no_population.tsv", checkIfExists: true)

    ch_input = ch_vcf
        .combine(ch_pop)
        .map { vcf, pop_file ->
            tuple([id: params.output_prefix], vcf, pop_file)
        }

    NJ_TREE(
        ch_input,
        params.output_prefix,
        params.window_size,
    )
}
