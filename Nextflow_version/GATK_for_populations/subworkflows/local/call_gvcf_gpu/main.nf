/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARABRICKS GPU PER-SAMPLE GVCF CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARABRICKS_HAPLOTYPECALLER } from '../../../modules/nf-core/parabricks/haplotypecaller/main'
include { PARABRICKS_INDEXGVCF       } from '../../../modules/nf-core/parabricks/indexgvcf/main'

workflow CALL_GVCF_GPU {

    take:
    ch_bam
    ch_bai
    ch_fasta
    intervals

    main:

    def ch_haplotypecaller_input = ch_bam
        .join(ch_bai)
        .map { meta, bam, bai -> [meta, bam, bai, intervals] }

    PARABRICKS_HAPLOTYPECALLER(ch_haplotypecaller_input, ch_fasta)
    PARABRICKS_INDEXGVCF(PARABRICKS_HAPLOTYPECALLER.out.gvcf)

    emit:
    gvcf = PARABRICKS_HAPLOTYPECALLER.out.gvcf
    tbi = PARABRICKS_INDEXGVCF.out.gvcf_index
}
