/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GATK4 CPU JOINT GENOTYPING BY CONTIG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_GENOTYPEGVCFS } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'

workflow JOINT_CALL_CPU {

    take:
    ch_genomicsdb
    ch_fasta
    ch_fai
    ch_dict

    main:

    def reference_meta = [id: 'reference']
    def ch_genotype_input = ch_genomicsdb.map { meta, genomicsdb ->
        [meta, genomicsdb, [], [], []]
    }
    def ch_empty_reference = channel.value([reference_meta, []])

    GATK4_GENOTYPEGVCFS(
        ch_genotype_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_empty_reference,
        ch_empty_reference
    )

    emit:
    vcf = GATK4_GENOTYPEGVCFS.out.vcf
    tbi = GATK4_GENOTYPEGVCFS.out.tbi
}
