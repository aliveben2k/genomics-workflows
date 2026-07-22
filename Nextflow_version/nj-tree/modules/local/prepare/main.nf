process PREPARE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/environment.yml"

    input:
    tuple val(meta), path(vcf), path(pop_file)
    val output_prefix
    val window_size

    output:
    tuple val(meta), path("${output_prefix}.*.tmp.txt.gz"), emit: chunks
    tuple val(meta), path("${output_prefix}.chunks.tsv"), emit: manifest

    script:
    def pop_arg = pop_file.name == 'no_population.tsv' ? '' : "--pop ${pop_file}"
    """
    vcf2table_large.py prepare \
        --input ${vcf} \
        --output ${output_prefix} \
        --window ${window_size} \
        ${pop_arg}
    """

    stub:
    """
    touch ${output_prefix}.1.tmp.txt.gz
    touch ${output_prefix}.chunks.tsv
    """
}
