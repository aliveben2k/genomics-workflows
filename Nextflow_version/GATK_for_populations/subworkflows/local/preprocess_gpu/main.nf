/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARABRICKS GPU READ ALIGNMENT AND GATK PREPROCESSING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARABRICKS_FQ2BAM     } from '../../../modules/nf-core/parabricks/fq2bam/main'
include { PARABRICKS_APPLYBQSR  } from '../../../modules/nf-core/parabricks/applybqsr/main'
include { FINALIZE_BQSR_BAM as FINALIZE_GPU_BQSR_BAM } from '../../../modules/local/finalize_bqsr_bam/main'

workflow PREPROCESS_GPU {

    take:
    ch_reads
    ch_bwa_index
    ch_fasta
    fasta
    known_sites
    intervals
    run_bqsr

    main:

    def reference_meta = [id: 'reference']
    def ch_intervals = channel.value([reference_meta, intervals])
    def ch_known_sites = channel.value([reference_meta, known_sites])

    PARABRICKS_FQ2BAM(
        ch_reads,
        ch_fasta,
        ch_bwa_index,
        ch_intervals,
        ch_known_sites,
        'bam'
    )

    def ch_final_bam = channel.empty()
    def ch_final_bai = channel.empty()
    def ch_bqsr_table = channel.empty()

    if (run_bqsr) {
        PARABRICKS_APPLYBQSR(
            PARABRICKS_FQ2BAM.out.bam,
            PARABRICKS_FQ2BAM.out.bai,
            PARABRICKS_FQ2BAM.out.bqsr_table,
            ch_intervals,
            ch_fasta
        )
        def ch_temporary_bqsr_bam = PARABRICKS_APPLYBQSR.out.bam.join(PARABRICKS_APPLYBQSR.out.bai)
        FINALIZE_GPU_BQSR_BAM(ch_temporary_bqsr_bam)

        ch_final_bam = FINALIZE_GPU_BQSR_BAM.out.bam
        ch_final_bai = FINALIZE_GPU_BQSR_BAM.out.bai
        ch_bqsr_table = PARABRICKS_FQ2BAM.out.bqsr_table
    } else {
        ch_final_bam = PARABRICKS_FQ2BAM.out.bam
        ch_final_bai = PARABRICKS_FQ2BAM.out.bai
    }

    emit:
    bam = ch_final_bam
    bai = ch_final_bai
    bqsr_table = ch_bqsr_table
    duplicate_metrics = PARABRICKS_FQ2BAM.out.duplicate_metrics
}
