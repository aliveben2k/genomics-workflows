/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GATK4 SAMPLE-SUBSET JOINT GENOTYPING BY CONTIG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CREATE_SUBSET_SAMPLE_ARGS               } from '../../../modules/local/create_subset_sample_args/main'
include { GATK4_SELECTVARIANTS_FROM_GENOMICSDB    } from '../../../modules/local/gatk4_selectvariants_from_genomicsdb/main'
include { GATK4_GENOTYPEGVCFS                     } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'

workflow JOINT_CALL_SUBSET_CPU {

    take:
    ch_genomicsdb
    ch_subset_sample_ids
    ch_fasta
    ch_fai
    ch_dict

    main:

    def reference_meta = [id: 'reference']
    def ch_empty_reference = channel.value([reference_meta, []])

    CREATE_SUBSET_SAMPLE_ARGS(ch_subset_sample_ids)

    GATK4_SELECTVARIANTS_FROM_GENOMICSDB(
        ch_genomicsdb,
        CREATE_SUBSET_SAMPLE_ARGS.out.sample_args.first(),
        ch_fasta,
        ch_fai,
        ch_dict
    )

    def ch_genotype_input = GATK4_SELECTVARIANTS_FROM_GENOMICSDB.out.gvcf
        .join(GATK4_SELECTVARIANTS_FROM_GENOMICSDB.out.tbi)
        .map { meta, gvcf, tbi -> [meta, gvcf, tbi, [], []] }

    GATK4_GENOTYPEGVCFS(
        ch_genotype_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_empty_reference,
        ch_empty_reference
    )

    emit:
    subset_gvcf = GATK4_SELECTVARIANTS_FROM_GENOMICSDB.out.gvcf
    subset_tbi  = GATK4_SELECTVARIANTS_FROM_GENOMICSDB.out.tbi
    vcf         = GATK4_GENOTYPEGVCFS.out.vcf
    tbi         = GATK4_GENOTYPEGVCFS.out.tbi
}
