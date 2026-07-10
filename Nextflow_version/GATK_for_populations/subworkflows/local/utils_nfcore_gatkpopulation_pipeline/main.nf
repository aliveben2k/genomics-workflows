//
// Subworkflow with functionality specific to the nf-core/gatkpopulation pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message
    start_stage       // integer: First stage to execute
    stop_stage        // integer: Last stage to execute

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(start_stage, stop_stage, input)

    //
    // Create channel from input file provided through params.input
    //

    ch_stage_input = channel.empty()
    if (input) {
        if (start_stage <= 2) {
            ch_stage_input = channel
                .fromPath(input, checkIfExists: true)
                .splitCsv(header: true, strip: true)
                .map { row ->
                    requireCsvColumns(row, ['sample', 'fastq_1', 'fastq_2'])
                    def sample_id = requireSimpleId(row.sample, 'sample')
                    def fastq_1 = requireInputFile(row.fastq_1, 'fastq_1', /^([\S\s]*\/)?[^\s\/]+\.f(ast)?q(\.trimmed)?\.gz$/)
                    def fastq_2 = row.fastq_2 ? requireInputFile(row.fastq_2, 'fastq_2', /^([\S\s]*\/)?[^\s\/]+\.f(ast)?q(\.trimmed)?\.gz$/) : null
                    [sample_id, [id: sample_id, single_end: !fastq_2], fastq_2 ? [fastq_1, fastq_2] : [fastq_1]]
                }
                .groupTuple()
                .map { samplesheet -> validateInputSamplesheet(samplesheet) }
                .map { meta, fastqs -> [meta, fastqs.flatten()] }
        } else if (start_stage == 3) {
            ch_stage_input = channel
                .fromPath(input, checkIfExists: true)
                .splitCsv(header: true, strip: true)
                .map { row ->
                    requireCsvColumns(row, ['sample', 'bam', 'bai'])
                    def sample_id = requireSimpleId(row.sample, 'sample')
                    def bam = requireInputFile(row.bam, 'bam', /^\S+\.bam$/)
                    def bai = requireInputFile(row.bai, 'bai', /^\S+\.(bam\.)?bai$/)
                    [sample_id, [id: sample_id], bam, bai]
                }
                .groupTuple()
                .map { id, metas, bams, bais ->
                    requireUniqueEntry(id, metas)
                    requireBamIndex(bams[0], bais[0])
                    [metas[0], bams[0], bais[0]]
                }
        } else if (start_stage == 4) {
            ch_stage_input = channel
                .fromPath(input, checkIfExists: true)
                .splitCsv(header: true, strip: true)
                .map { row ->
                    requireCsvColumns(row, ['sample', 'gvcf', 'tbi'])
                    def sample_id = requireSimpleId(row.sample, 'sample')
                    def gvcf = requireInputFile(row.gvcf, 'gvcf', /^\S+\.g\.vcf\.gz$/)
                    def tbi = requireInputFile(row.tbi, 'tbi', /^\S+\.g\.vcf\.gz\.tbi$/)
                    [sample_id, [id: sample_id], gvcf, tbi]
                }
                .groupTuple()
                .map { id, metas, gvcfs, tbis ->
                    requireUniqueEntry(id, metas)
                    requireTabixIndex(gvcfs[0], tbis[0])
                    [metas[0], gvcfs[0], tbis[0]]
                }
        } else if (start_stage == 5) {
            ch_stage_input = channel
                .fromPath(input, checkIfExists: true)
                .splitCsv(header: true, strip: true)
                .map { row ->
                    requireCsvColumns(row, ['contig', 'genomicsdb'])
                    def contig = requireSimpleId(row.contig, 'contig')
                    [contig, [id: contig], file(row.genomicsdb, checkIfExists: true)]
                }
                .groupTuple()
                .map { id, metas, databases ->
                    requireUniqueEntry(id, metas)
                    [metas[0] + [contig: id], databases[0]]
                }
        } else {
            ch_stage_input = channel
                .fromPath(input, checkIfExists: true)
                .splitCsv(header: true, strip: true)
                .map { row ->
                    requireCsvColumns(row, ['contig', 'vcf', 'tbi'])
                    def contig = requireSimpleId(row.contig, 'contig')
                    def vcf = requireInputFile(row.vcf, 'vcf', /^\S+\.vcf\.gz$/)
                    def tbi = requireInputFile(row.tbi, 'tbi', /^\S+\.vcf\.gz\.tbi$/)
                    [contig, [id: contig], vcf, tbi]
                }
                .groupTuple()
                .map { id, metas, vcfs, tbis ->
                    requireUniqueEntry(id, metas)
                    requireTabixIndex(vcfs[0], tbis[0])
                    [metas[0] + [contig: id], vcfs[0], tbis[0]]
                }
        }
    }

    emit:
    stage_input = ch_stage_input
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters(start_stage, stop_stage, input) {
    genomeExistsError()
    if (start_stage < 0 || start_stage > 6) {
        error '--start_stage must be between 0 and 6.'
    }
    if (stop_stage < start_stage || stop_stage > 6) {
        error '--stop_stage must be between --start_stage and 6.'
    }
    if (!input && !(start_stage == 0 && stop_stage == 0)) {
        error "--input is required for --start_stage ${start_stage}."
    }
}

def requireUniqueEntry(id, metas) {
    if (metas.size() != 1) {
        error "Entry IDs must be unique in the stage input manifest: ${id}"
    }
}

def requireCsvColumns(row, expected_columns) {
    def observed_columns = row.keySet() as Set
    def expected_column_set = expected_columns as Set
    if (observed_columns != expected_column_set) {
        error "Input CSV must contain exactly these columns: ${expected_columns.join(', ')}"
    }
}

def requireSimpleId(value, column_name) {
    def id = value?.toString()?.trim()
    if (!id || id =~ /\s/) {
        error "Column '${column_name}' must contain a non-empty ID without whitespace."
    }
    id
}

def requireInputFile(value, column_name, pattern) {
    def path = value?.toString()?.trim()
    if (!path) {
        error "Column '${column_name}' must contain a file path."
    }
    if (!(path ==~ pattern)) {
        error "Column '${column_name}' has an invalid file name: ${path}"
    }
    file(path, checkIfExists: true)
}

def requireBamIndex(bam, bai) {
    def valid_names = ["${bam.name}.bai", "${bam.baseName}.bai"]
    if (!valid_names.contains(bai.name)) {
        error "BAM index ${bai.name} does not match ${bam.name}."
    }
}

def requireTabixIndex(vcf, tbi) {
    if (tbi.name != "${vcf.name}.tbi") {
        error "Tabix index ${tbi.name} does not match ${vcf.name}."
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Lane merging will be added as an explicit process; until then each sample is one row.
    if (metas.size() != 1) {
        error("Please check input samplesheet -> Sample IDs must be unique (one row per sample): ${metas[0].id}")
    }

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "fastp (Chen et al. 2018),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Chen S, Zhou Y, Chen Y, Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. doi: 10.1093/bioinformatics/bty560.</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
