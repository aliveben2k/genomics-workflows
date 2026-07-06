/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GATK4 CPU PER-SAMPLE GVCF CALLING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow CALL_GVCF_CPU {

    take:
    ch_bam
    ch_bai
    ch_fasta
    ch_fasta_fai
    ch_fasta_dict
    interval

    main:

    def ch_haplotypecaller_input = ch_bam
        .join(ch_bai)
        .map { meta, bam, bai -> [meta, bam, bai, interval, []] }
    def reference_meta = [id: 'reference']
    def ch_empty_reference = channel.value([reference_meta, []])

    GATK4_HAPLOTYPECALLER(
        ch_haplotypecaller_input,
        ch_fasta,
        ch_fasta_fai,
        ch_fasta_dict,
        ch_empty_reference,
        ch_empty_reference
    )

    emit:
    gvcf = GATK4_HAPLOTYPECALLER.out.vcf
    tbi = GATK4_HAPLOTYPECALLER.out.tbi
}
