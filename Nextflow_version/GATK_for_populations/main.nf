#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/gatkpopulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/gatkpopulation
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATKPOPULATION  } from './workflows/gatkpopulation'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_gatkpopulation_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_gatkpopulation_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_gatkpopulation_pipeline'

def runSerialFromSession(session_id) {
    def encoded = new BigInteger(
        session_id.toString().replace('-', ''),
        16
    ).toString(36).toUpperCase()
    def first = ('A'..'Z')[Math.floorMod(encoded.hashCode(), 26)]
    "${first}${encoded.take(3)}"
}

def asBooleanParam(value) {
    if (value instanceof Boolean) {
        return value
    }
    return value?.toString()?.toBoolean()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_GATKPOPULATION {

    take:
    stage_input // channel: stage-specific records read from --input
    start_stage
    stop_stage
    fasta
    gather_vcfs
    sequencing_type

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GATKPOPULATION (
        stage_input,
        start_stage,
        stop_stage,
        params.trimmer,
        params.adapter_fasta,
        fasta,
        params.fasta_fai,
        params.fasta_dict,
        params.bwa_index,
        params.known_sites,
        params.intervals,
        asBooleanParam(params.gpu_fallback),
        sequencing_type,
        params.all_sites,
        params.subset_samples,
        gather_vcfs,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    processed_reads = GATKPOPULATION.out.processed_reads // channel: [ val(meta), path(reads) ]
    preprocessed_bam = GATKPOPULATION.out.preprocessed_bam
    preprocessed_bai = GATKPOPULATION.out.preprocessed_bai
    bqsr_table       = GATKPOPULATION.out.bqsr_table
    bam_stats        = GATKPOPULATION.out.bam_stats
    sample_gvcf      = GATKPOPULATION.out.sample_gvcf
    sample_gvcf_tbi  = GATKPOPULATION.out.sample_gvcf_tbi
    genomicsdb       = GATKPOPULATION.out.genomicsdb
    raw_vcf          = GATKPOPULATION.out.raw_vcf
    raw_vcf_tbi      = GATKPOPULATION.out.raw_vcf_tbi
    filtered_vcf     = GATKPOPULATION.out.filtered_vcf
    filtered_vcf_tbi = GATKPOPULATION.out.filtered_vcf_tbi
    multiqc_report  = GATKPOPULATION.out.multiqc_report  // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    // TODO nf-core: Remove this line if you don't need a FASTA file
    //   This is an example of how to use getGenomeAttribute() to fetch parameters
    //   from igenomes.config using `--genome`
    def fasta = params.fasta ?: getGenomeAttribute('fasta')
    def start_stage = params.start_stage as Integer
    def stop_stage = params.stop_stage == null ? 6 : params.stop_stage as Integer
    def gather_vcfs = params.gather_vcfs == null ? null : params.gather_vcfs as Integer
    def filter_vcf = params.filter_vcf ?: 'biallele'
    def sequencing_type = (params.sequencing_type ?: 'wgs').toString().toLowerCase()

    def sequencing_type_aliases = ['wgs', 'wes', 'rad', 'amp'].findAll { type -> asBooleanParam(params[type]) }
    if (sequencing_type_aliases.size() > 1) {
        error "Use only one sequencing-type flag: ${sequencing_type_aliases.collect { '--' + it }.join(', ')}"
    }
    if (!['wgs', 'wes', 'rad', 'amp'].contains(sequencing_type)) {
        error "--sequencing_type must be 'wgs', 'wes', 'rad', or 'amp'."
    }
    if (!sequencing_type_aliases.isEmpty()) {
        def alias_type = sequencing_type_aliases[0]
        if (sequencing_type != 'wgs' && sequencing_type != alias_type) {
            error "--${alias_type} conflicts with --sequencing_type ${sequencing_type}."
        }
        sequencing_type = alias_type
    }
    if (!['biallele', 'monobi'].contains(filter_vcf)) {
        error "--filter_vcf must be 'biallele' or 'monobi'."
    }
    if (filter_vcf == 'monobi' && params.bcftools_view_args == null && !params.all_sites && start_stage <= 5) {
        log.warn '--filter_vcf monobi cannot retain monomorphic sites unless joint genotyping uses --all_sites.'
    }
    if (gather_vcfs != null && ![5, 6].contains(gather_vcfs)) {
        error '--gather_vcfs must be 5 or 6.'
    }
    if (gather_vcfs != null && stop_stage < gather_vcfs) {
        error "--gather_vcfs ${gather_vcfs} requires --stop_stage ${gather_vcfs} or later."
    }
    if (gather_vcfs == 5 && start_stage > 5) {
        error '--gather_vcfs 5 requires --start_stage 5 or earlier. To filter an existing all.raw.vcf.gz, provide it as the stage-6 input and omit --gather_vcfs 5.'
    }
    def run_serial = params.run_serial ?: runSerialFromSession(workflow.sessionId)

    log.info "Pipeline output serial: ${run_serial}"
    log.info "Pipeline stage range: ${start_stage}-${stop_stage}"
    log.info "Sequencing type: ${sequencing_type}; " +
        "duplicate marking: ${['wgs', 'wes'].contains(sequencing_type) ? 'enabled' : 'disabled'}"

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        start_stage,
        stop_stage
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_GATKPOPULATION (
        PIPELINE_INITIALISATION.out.stage_input,
        start_stage,
        stop_stage,
        fasta,
        gather_vcfs,
        sequencing_type
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_GATKPOPULATION.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
