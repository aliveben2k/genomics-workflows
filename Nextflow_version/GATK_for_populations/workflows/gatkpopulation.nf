/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { TRIMMOMATIC            } from '../modules/nf-core/trimmomatic/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_STATS         } from '../modules/nf-core/samtools/stats/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_GENOMICSDBIMPORT } from '../modules/nf-core/gatk4/genomicsdbimport/main'
include { BCFTOOLS_VIEW          } from '../modules/nf-core/bcftools/view/main'
include { GATK4_GATHERVCFS as GATHER_RAW_VCFS      } from '../modules/local/gatk4_gathervcfs/main'
include { GATK4_GATHERVCFS as GATHER_FILTERED_VCFS } from '../modules/local/gatk4_gathervcfs/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { CREATE_GVCF_SAMPLE_MAP } from '../modules/local/create_gvcf_sample_map/main'
include { PREPROCESS_CPU         } from '../subworkflows/local/preprocess_cpu/main'
include { PREPROCESS_GPU         } from '../subworkflows/local/preprocess_gpu/main'
include { CALL_GVCF_CPU          } from '../subworkflows/local/call_gvcf_cpu/main'
include { CALL_GVCF_GPU          } from '../subworkflows/local/call_gvcf_gpu/main'
include { JOINT_CALL_CPU         } from '../subworkflows/local/joint_call_cpu/main'
include { JOINT_CALL_SUBSET_CPU  } from '../subworkflows/local/joint_call_subset_cpu/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_gatkpopulation_pipeline'

def prepareGatherInput(records, contig_order, source_label) {
    if (records.isEmpty()) {
        error "No ${source_label} VCFs were available for GatherVcfs."
    }
    def contig_rank = contig_order.withIndex().collectEntries { contig, index -> [(contig): index] }
    def record_contigs = records.collect { meta, _vcf, _tbi -> meta.contig ?: meta.id }
    def unknown_contigs = record_contigs.findAll { !contig_rank.containsKey(it) }.unique()
    if (!unknown_contigs.isEmpty()) {
        error "Cannot order ${source_label} VCFs because these contigs are absent from the reference FASTA index: ${unknown_contigs.sort().join(', ')}"
    }
    def duplicate_contigs = record_contigs.countBy { it }.findAll { _contig, count -> count > 1 }.keySet()
    if (!duplicate_contigs.isEmpty()) {
        error "Cannot gather duplicate ${source_label} contigs: ${duplicate_contigs.sort().join(', ')}"
    }
    def sorted_records = records.sort { left, right ->
        def left_contig = left[0].contig ?: left[0].id
        def right_contig = right[0].contig ?: right[0].id
        contig_rank[left_contig] <=> contig_rank[right_contig]
    }
    [
        [id: 'all', contig: 'all'],
        sorted_records.collect { _meta, vcf, _tbi -> vcf },
        sorted_records.collect { _meta, _vcf, tbi -> tbi }
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GATKPOPULATION {

    take:
    ch_stage_input // channel: records matching --start_stage
    start_stage
    stop_stage
    trimmer
    adapter_fasta
    fasta
    fasta_fai
    fasta_dict
    bwa_index
    known_sites
    intervals
    gpu_fallback
    sequencing_type
    all_sites
    subset_samples
    gather_vcfs
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def mark_duplicates = ['wgs', 'wes'].contains(sequencing_type)
    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def ch_raw_reads = channel.empty()
    def ch_trimmed_reads = channel.empty()
    def reference_meta = [id: 'reference']
    def fasta_file = []
    def ch_fasta = channel.empty()
    def ch_reference_fai = channel.empty()
    def ch_reference_dict = channel.empty()
    def fai_file = []
    def dict_file = []
    def ch_bwa_index = channel.empty()
    def ch_preprocessed_bam = channel.empty()
    def ch_preprocessed_bai = channel.empty()
    def ch_bqsr_table = channel.empty()
    def ch_duplicate_metrics = channel.empty()
    def ch_bam_stats = channel.empty()
    def ch_sample_gvcf = channel.empty()
    def ch_sample_gvcf_tbi = channel.empty()
    def ch_genomicsdb = channel.empty()
    def ch_raw_vcf = channel.empty()
    def ch_raw_vcf_tbi = channel.empty()
    def ch_filtered_vcf = channel.empty()
    def ch_filtered_vcf_tbi = channel.empty()
    def ch_subset_sample_ids = channel.value([])
    def ch_reference_contig_order = channel.empty()

    if (subset_samples && start_stage <= 5 && stop_stage >= 5) {
        ch_subset_sample_ids = channel
            .fromPath(subset_samples, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { row ->
                if ((row.keySet() as Set) != (['sample'] as Set)) {
                    error '--subset_samples must contain exactly one CSV column named sample.'
                }
                def sample_id = row.sample?.toString()?.trim()
                if (!sample_id || sample_id =~ /\s/) {
                    error "Invalid sample ID '${sample_id ?: ''}' in --subset_samples; IDs must be non-empty and contain no whitespace."
                }
                sample_id
            }
            .collect()
            .map { sample_ids ->
                if (sample_ids.isEmpty()) {
                    error '--subset_samples must contain at least one sample.'
                }
                def duplicate_ids = sample_ids.countBy { it }.findAll { _id, count -> count > 1 }.keySet()
                if (!duplicate_ids.isEmpty()) {
                    error "Duplicate sample IDs in --subset_samples: ${duplicate_ids.sort().join(', ')}"
                }
                sample_ids
            }
    }

    if (start_stage <= 1) {
        ch_raw_reads = ch_stage_input
    } else if (start_stage == 2) {
        ch_trimmed_reads = ch_stage_input
    } else if (start_stage == 3) {
        ch_preprocessed_bam = ch_stage_input.map { meta, bam, _bai -> [meta, bam] }
        ch_preprocessed_bai = ch_stage_input.map { meta, _bam, bai -> [meta, bai] }
    } else if (start_stage == 4) {
        ch_sample_gvcf = ch_stage_input.map { meta, gvcf, _tbi -> [meta, gvcf] }
        ch_sample_gvcf_tbi = ch_stage_input.map { meta, _gvcf, tbi -> [meta, tbi] }
    } else if (start_stage == 5) {
        ch_genomicsdb = ch_stage_input
    } else if (start_stage == 6) {
        ch_raw_vcf = ch_stage_input.map { meta, vcf, _tbi -> [meta, vcf] }
        ch_raw_vcf_tbi = ch_stage_input.map { meta, _vcf, tbi -> [meta, tbi] }
    }

    def needs_reference = start_stage == 0 || (start_stage <= 5 && stop_stage >= 2)
    if (needs_reference) {
        if (!fasta) {
            error 'A reference FASTA is required for the selected stage range.'
        }
        fasta_file = file(fasta, checkIfExists: true)
        if (fasta_file.name.endsWith('.gz')) {
            error 'The pipeline requires an uncompressed reference FASTA.'
        }
        ch_fasta = channel.value([reference_meta, fasta_file])

        if (start_stage == 0) {
            SAMTOOLS_FAIDX(channel.value([reference_meta, fasta_file, []]), false)
            GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
            ch_reference_fai = SAMTOOLS_FAIDX.out.fai
            ch_reference_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
            if (bwa_index) {
                ch_bwa_index = channel.value([reference_meta, file(bwa_index, checkIfExists: true)])
            } else {
                BWA_INDEX(ch_fasta)
                ch_bwa_index = BWA_INDEX.out.index.first()
            }
        } else {
            def inferred_fai = fasta_fai ?: "${fasta_file}.fai"
            def inferred_dict = fasta_dict ?: "${fasta_file.parent}/${fasta_file.baseName}.dict"
            ch_reference_fai = channel.value([reference_meta, file(inferred_fai, checkIfExists: true)])
            ch_reference_dict = channel.value([reference_meta, file(inferred_dict, checkIfExists: true)])
            if (start_stage <= 2 && stop_stage >= 2) {
                if (!bwa_index) {
                    error '--bwa_index is required when starting at stage 1 or 2 and running stage 2.'
                }
                ch_bwa_index = channel.value([reference_meta, file(bwa_index, checkIfExists: true)])
            }
        }
        fai_file = ch_reference_fai.map { _meta, fai -> fai }.first()
        dict_file = ch_reference_dict.map { _meta, dict -> dict }.first()
        if (gather_vcfs != null) {
            ch_reference_contig_order = ch_reference_fai
                .map { _meta, fai -> fai }
                .splitText()
                .map { line -> line.tokenize('\t')[0] }
                .collect()
        }
    }
    if (start_stage == 6 && gather_vcfs == 6) {
        if (!fasta_fai) {
            error '--fasta_fai is required when starting at stage 6 with --gather_vcfs 6.'
        }
        ch_reference_contig_order = channel
            .fromPath(fasta_fai, checkIfExists: true)
            .splitText()
            .map { line -> line.tokenize('\t')[0] }
            .collect()
    }

    def uses_bam_or_gvcf_calling = start_stage <= 3 && stop_stage >= 2
    def known_site_paths = !uses_bam_or_gvcf_calling || !known_sites
        ? []
        : known_sites instanceof List
            ? known_sites
            : known_sites.toString().split(',').collect { it.trim() }.findAll { it }
    def known_site_files = known_site_paths.collect { path ->
        file(path, checkIfExists: true)
    }
    def known_site_tbi_files = known_site_paths.collect { path ->
        file("${path}.tbi", checkIfExists: true)
    }
    def run_bqsr = !known_site_files.isEmpty()
    if (start_stage <= 2 && stop_stage >= 2 && !run_bqsr) {
        log.warn 'No --known_sites were supplied; BQSR will be skipped.'
    }
    def interval_file = uses_bam_or_gvcf_calling && intervals ? file(intervals, checkIfExists: true) : []
    def interval_files = interval_file ? [interval_file] : []

    //
    // STAGE 1: Adapter/quality trimming and FastQC
    //
    if (start_stage <= 1 && stop_stage >= 1) {
        if (trimmer == 'fastp') {
            FASTP(
                ch_raw_reads.map { meta, reads ->
                    [meta, reads, adapter_fasta ? file(adapter_fasta, checkIfExists: true) : []]
                },
                false,
                false,
                false
            )
            ch_trimmed_reads = FASTP.out.reads
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { _meta, report -> report })
        } else if (trimmer == 'trimmomatic') {
            TRIMMOMATIC(ch_raw_reads)
            ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
            ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.out_log.map { _meta, report -> report })
        } else {
            error "Invalid trimmer '${trimmer}'. Choose 'fastp' or 'trimmomatic'."
        }
        FASTQC(ch_trimmed_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map { _meta, report -> report })
    }

    //
    // STAGE 2: BAM preprocessing and statistics
    //
    if (start_stage <= 2 && stop_stage >= 2) {
        if (workflow.profile.tokenize(',').contains('gpu')) {
            PREPROCESS_GPU(
                ch_trimmed_reads,
                ch_bwa_index,
                ch_fasta,
                fasta_file,
                known_site_files,
                interval_files,
                run_bqsr
            )

            if (gpu_fallback) {
                def ch_gpu_results = PREPROCESS_GPU.out.bam.join(PREPROCESS_GPU.out.bai)
                def ch_gpu_success_ids = ch_gpu_results
                    .map { meta, _bam, _bai -> meta.id }
                    .collect()
                    .map { ids -> [ids] }
                def ch_cpu_fallback_reads = ch_trimmed_reads
                    .combine(ch_gpu_success_ids)
                    .filter { meta, _reads, successful_ids -> !successful_ids.contains(meta.id) }
                    .map { meta, reads, _successful_ids ->
                        log.warn "GPU preprocessing failed or was unavailable for ${meta.id}; falling back to CPU."
                        [meta, reads]
                    }

                PREPROCESS_CPU(
                    ch_cpu_fallback_reads,
                    ch_bwa_index,
                    ch_fasta,
                    fasta_file,
                    fai_file,
                    dict_file,
                    known_site_files,
                    known_site_tbi_files,
                    interval_file,
                    run_bqsr,
                    mark_duplicates
                )

                def ch_successful_gpu_bam = ch_gpu_results.map { meta, bam, _bai -> [meta, bam] }
                def ch_successful_gpu_bai = ch_gpu_results.map { meta, _bam, bai -> [meta, bai] }
                def ch_successful_gpu_tables = PREPROCESS_GPU.out.bqsr_table
                    .combine(ch_gpu_success_ids)
                    .filter { meta, _table, successful_ids -> successful_ids.contains(meta.id) }
                    .map { meta, table, _successful_ids -> [meta, table] }
                def ch_successful_gpu_metrics = PREPROCESS_GPU.out.duplicate_metrics
                    .combine(ch_gpu_success_ids)
                    .filter { meta, _metrics, successful_ids -> successful_ids.contains(meta.id) }
                    .map { meta, metrics, _successful_ids -> [meta, metrics] }
                ch_preprocessed_bam = ch_successful_gpu_bam.mix(PREPROCESS_CPU.out.bam)
                ch_preprocessed_bai = ch_successful_gpu_bai.mix(PREPROCESS_CPU.out.bai)
                ch_bqsr_table = ch_successful_gpu_tables.mix(PREPROCESS_CPU.out.bqsr_table)
                ch_duplicate_metrics = ch_successful_gpu_metrics.mix(PREPROCESS_CPU.out.duplicate_metrics)
            } else {
                ch_preprocessed_bam = PREPROCESS_GPU.out.bam
                ch_preprocessed_bai = PREPROCESS_GPU.out.bai
                ch_bqsr_table = PREPROCESS_GPU.out.bqsr_table
                ch_duplicate_metrics = PREPROCESS_GPU.out.duplicate_metrics
            }
        } else {
            PREPROCESS_CPU(
                ch_trimmed_reads,
                ch_bwa_index,
                ch_fasta,
                fasta_file,
                fai_file,
                dict_file,
                known_site_files,
                known_site_tbi_files,
                interval_file,
                run_bqsr,
                mark_duplicates
            )
            ch_preprocessed_bam = PREPROCESS_CPU.out.bam
            ch_preprocessed_bai = PREPROCESS_CPU.out.bai
            ch_bqsr_table = PREPROCESS_CPU.out.bqsr_table
            ch_duplicate_metrics = PREPROCESS_CPU.out.duplicate_metrics
        }
        SAMTOOLS_STATS(
            ch_preprocessed_bam.join(ch_preprocessed_bai),
            ch_reference_fai.map { meta, fai -> [meta, fasta_file, fai] }
        )
        ch_bam_stats = SAMTOOLS_STATS.out.stats
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_duplicate_metrics.map { _meta, report -> report })
            .mix(ch_bam_stats.map { _meta, report -> report })
    }

    //
    // STAGE 3: Per-sample GVCF calling
    //
    if (start_stage <= 3 && stop_stage >= 3 && workflow.profile.tokenize(',').contains('gpu')) {
        CALL_GVCF_GPU(
            ch_preprocessed_bam,
            ch_preprocessed_bai,
            ch_fasta,
            interval_files
        )

        if (gpu_fallback) {
            def ch_gpu_gvcf_results = CALL_GVCF_GPU.out.gvcf.join(CALL_GVCF_GPU.out.tbi)
            def ch_gpu_gvcf_success_ids = ch_gpu_gvcf_results
                .map { meta, _gvcf, _tbi -> meta.id }
                .collect()
                .map { ids -> [ids] }
            def ch_cpu_fallback_bam = ch_preprocessed_bam
                .combine(ch_gpu_gvcf_success_ids)
                .filter { meta, _bam, successful_ids -> !successful_ids.contains(meta.id) }
                .map { meta, bam, _successful_ids ->
                    log.warn "GPU GVCF calling failed or was unavailable for ${meta.id}; falling back to CPU."
                    [meta, bam]
                }
            def ch_cpu_fallback_bai = ch_preprocessed_bai
                .combine(ch_gpu_gvcf_success_ids)
                .filter { meta, _bai, successful_ids -> !successful_ids.contains(meta.id) }
                .map { meta, bai, _successful_ids -> [meta, bai] }

            CALL_GVCF_CPU(
                ch_cpu_fallback_bam,
                ch_cpu_fallback_bai,
                ch_fasta,
                ch_reference_fai,
                ch_reference_dict,
                interval_file
            )

            def ch_successful_gpu_gvcf = ch_gpu_gvcf_results.map { meta, gvcf, _tbi -> [meta, gvcf] }
            def ch_successful_gpu_tbi = ch_gpu_gvcf_results.map { meta, _gvcf, tbi -> [meta, tbi] }
            ch_sample_gvcf = ch_successful_gpu_gvcf.mix(CALL_GVCF_CPU.out.gvcf)
            ch_sample_gvcf_tbi = ch_successful_gpu_tbi.mix(CALL_GVCF_CPU.out.tbi)
        } else {
            ch_sample_gvcf = CALL_GVCF_GPU.out.gvcf
            ch_sample_gvcf_tbi = CALL_GVCF_GPU.out.tbi
        }
    } else if (start_stage <= 3 && stop_stage >= 3) {
        CALL_GVCF_CPU(
            ch_preprocessed_bam,
            ch_preprocessed_bai,
            ch_fasta,
            ch_reference_fai,
            ch_reference_dict,
            interval_file
        )
        ch_sample_gvcf = CALL_GVCF_CPU.out.gvcf
        ch_sample_gvcf_tbi = CALL_GVCF_CPU.out.tbi
    }

    //
    // STAGE 4: One all-sample GenomicsDB per contig
    //
    if (start_stage <= 4 && stop_stage >= 4) {
        def ch_cohort_gvcfs = ch_sample_gvcf
            .join(ch_sample_gvcf_tbi)
            .collect()
            .map { records ->
                if (records.isEmpty()) {
                    error 'No indexed sample GVCFs were available for GenomicsDBImport.'
                }
                def sorted_records = records.sort { left, right -> left[0].id <=> right[0].id }
                [
                    sorted_records.collect { meta, _gvcf, _tbi -> meta.id },
                    sorted_records.collect { _meta, gvcf, _tbi -> gvcf },
                    sorted_records.collect { _meta, _gvcf, tbi -> tbi }
                ]
            }

        CREATE_GVCF_SAMPLE_MAP(
            ch_cohort_gvcfs.map { sample_ids, gvcfs, _tbis -> [sample_ids, gvcfs] }
        )

        def ch_cohort_with_map = ch_cohort_gvcfs.combine(CREATE_GVCF_SAMPLE_MAP.out.sample_map)
        def ch_reference_contigs = ch_reference_fai
            .map { _meta, fai -> fai }
            .splitText()
            .map { line -> line.tokenize('\t')[0] }
        def ch_genomicsdb_input = ch_reference_contigs
            .combine(ch_cohort_with_map)
            .map { contig, _sample_ids, gvcfs, tbis, sample_map ->
                def safe_contig = contig.replaceAll(/[^A-Za-z0-9._-]/, '_')
                def workspace_id = safe_contig == contig
                    ? contig
                    : "${safe_contig}_${Integer.toUnsignedString(contig.hashCode(), 36)}"
                [[id: workspace_id, contig: contig], [sample_map] + gvcfs, tbis, [], contig, []]
            }

        GATK4_GENOMICSDBIMPORT(ch_genomicsdb_input, false, false, true)
        ch_genomicsdb = GATK4_GENOMICSDBIMPORT.out.genomicsdb
    }

    //
    // STAGE 5: GATK4 joint genotyping by contig
    //
    if (start_stage <= 5 && stop_stage >= 5) {
        if (subset_samples) {
            JOINT_CALL_SUBSET_CPU(
                ch_genomicsdb,
                ch_subset_sample_ids,
                ch_fasta,
                ch_reference_fai,
                ch_reference_dict
            )
            ch_raw_vcf = JOINT_CALL_SUBSET_CPU.out.vcf
            ch_raw_vcf_tbi = JOINT_CALL_SUBSET_CPU.out.tbi
        } else {
            JOINT_CALL_CPU(ch_genomicsdb, ch_fasta, ch_reference_fai, ch_reference_dict)
            ch_raw_vcf = JOINT_CALL_CPU.out.vcf
            ch_raw_vcf_tbi = JOINT_CALL_CPU.out.tbi
        }

        if (gather_vcfs == 5) {
            def ch_raw_gather_input = ch_raw_vcf
                .join(ch_raw_vcf_tbi)
                .collect()
                .map { records -> [records] }
                .combine(ch_reference_contig_order.map { contig_order -> [contig_order] })
                .map { records, contig_order ->
                    prepareGatherInput(records, contig_order, 'raw')
                }
            GATHER_RAW_VCFS(ch_raw_gather_input)
            ch_raw_vcf = GATHER_RAW_VCFS.out.vcf
            ch_raw_vcf_tbi = GATHER_RAW_VCFS.out.tbi
        }
    }

    //
    // STAGE 6: Filter each per-contig multi-sample VCF
    //
    if (start_stage <= 6 && stop_stage >= 6) {
        BCFTOOLS_VIEW(ch_raw_vcf.join(ch_raw_vcf_tbi), [], [], [])
        ch_filtered_vcf = BCFTOOLS_VIEW.out.vcf
        ch_filtered_vcf_tbi = BCFTOOLS_VIEW.out.index

        if (gather_vcfs == 6) {
            def ch_filtered_gather_input = ch_filtered_vcf
                .join(ch_filtered_vcf_tbi)
                .collect()
                .map { records -> [records] }
                .combine(ch_reference_contig_order.map { contig_order -> [contig_order] })
                .map { records, contig_order ->
                    prepareGatherInput(records, contig_order, 'filtered')
                }
            GATHER_FILTERED_VCFS(ch_filtered_gather_input)
            ch_filtered_vcf = GATHER_FILTERED_VCFS.out.vcf
            ch_filtered_vcf_tbi = GATHER_FILTERED_VCFS.out.tbi
        }
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'gatkpopulation_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'gatkpopulation'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    processed_reads = ch_trimmed_reads                                             // channel: [ val(meta), path(reads) ]
    preprocessed_bam = ch_preprocessed_bam                                         // channel: [ val(meta), path(bam) ]
    preprocessed_bai = ch_preprocessed_bai                                         // channel: [ val(meta), path(bai) ]
    bqsr_table = ch_bqsr_table                                                     // channel: [ val(meta), path(table) ]
    bam_stats = ch_bam_stats                                                        // channel: [ val(meta), path(stats) ]
    sample_gvcf = ch_sample_gvcf                                                    // channel: [ val(meta), path(g.vcf.gz) ]
    sample_gvcf_tbi = ch_sample_gvcf_tbi                                            // channel: [ val(meta), path(g.vcf.gz.tbi) ]
    genomicsdb = ch_genomicsdb                                                      // channel: [ val(meta), path(genomicsdb) ]
    raw_vcf = ch_raw_vcf                                                            // channel: [ val(meta), path(raw.vcf.gz) ]
    raw_vcf_tbi = ch_raw_vcf_tbi                                                    // channel: [ val(meta), path(raw.vcf.gz.tbi) ]
    filtered_vcf = ch_filtered_vcf                                                  // channel: [ val(meta), path(filtered.vcf.gz) ]
    filtered_vcf_tbi = ch_filtered_vcf_tbi                                          // channel: [ val(meta), path(filtered.vcf.gz.tbi) ]
    multiqc_report  = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions        = ch_versions                                                   // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
