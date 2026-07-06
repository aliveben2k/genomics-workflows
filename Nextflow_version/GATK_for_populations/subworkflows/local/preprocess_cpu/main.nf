/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPU READ ALIGNMENT AND GATK PREPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_MEM                  } from '../../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX           } from '../../../modules/local/samtools_index/main'
include { GATK4_MARKDUPLICATES     } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR   } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR          } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { FINALIZE_BQSR_BAM as FINALIZE_CPU_BQSR_BAM } from '../../../modules/local/finalize_bqsr_bam/main'

workflow PREPROCESS_CPU {

    take:
    ch_reads
    ch_bwa_index
    ch_fasta
    fasta
    fai
    dict
    known_sites
    known_sites_tbi
    interval
    run_bqsr
    mark_duplicates

    main:

    def reference_meta = [id: 'reference']
    def ch_reference_fasta = channel.value([reference_meta, fasta])
    def ch_reference_fai = channel.value([reference_meta, fai])
    def ch_reference_dict = channel.value([reference_meta, dict])
    def ch_known_sites = channel.value([reference_meta, known_sites])
    def ch_known_sites_tbi = channel.value([reference_meta, known_sites_tbi])

    BWA_MEM(ch_reads, ch_bwa_index, ch_fasta, true)

    def ch_pre_bqsr_bam = channel.empty()
    def ch_duplicate_metrics = channel.empty()
    if (mark_duplicates) {
        GATK4_MARKDUPLICATES(BWA_MEM.out.bam, fasta, fai)
        ch_pre_bqsr_bam = GATK4_MARKDUPLICATES.out.bam.join(GATK4_MARKDUPLICATES.out.bai)
        ch_duplicate_metrics = GATK4_MARKDUPLICATES.out.metrics
    } else {
        SAMTOOLS_INDEX(BWA_MEM.out.bam)
        ch_pre_bqsr_bam = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)
    }

    def ch_final_bam = channel.empty()
    def ch_final_bai = channel.empty()
    def ch_bqsr_table = channel.empty()

    if (run_bqsr) {
        def ch_baserecalibrator_input = ch_pre_bqsr_bam.map { meta, bam, bai ->
            [meta, bam, bai, interval]
        }
        GATK4_BASERECALIBRATOR(
            ch_baserecalibrator_input,
            ch_reference_fasta,
            ch_reference_fai,
            ch_reference_dict,
            ch_known_sites,
            ch_known_sites_tbi
        )

        def ch_applybqsr_input = ch_pre_bqsr_bam
            .join(GATK4_BASERECALIBRATOR.out.table)
            .map { meta, bam, bai, table -> [meta, bam, bai, table, interval] }
        GATK4_APPLYBQSR(ch_applybqsr_input, fasta, fai, dict)

        def ch_temporary_bqsr_bam = GATK4_APPLYBQSR.out.bam.join(GATK4_APPLYBQSR.out.bai)
        FINALIZE_CPU_BQSR_BAM(ch_temporary_bqsr_bam)

        ch_final_bam = FINALIZE_CPU_BQSR_BAM.out.bam
        ch_final_bai = FINALIZE_CPU_BQSR_BAM.out.bai
        ch_bqsr_table = GATK4_BASERECALIBRATOR.out.table
    } else {
        ch_final_bam = ch_pre_bqsr_bam.map { meta, bam, _bai -> [meta, bam] }
        ch_final_bai = ch_pre_bqsr_bam.map { meta, _bam, bai -> [meta, bai] }
    }

    emit:
    bam = ch_final_bam
    bai = ch_final_bai
    bqsr_table = ch_bqsr_table
    duplicate_metrics = ch_duplicate_metrics
}
