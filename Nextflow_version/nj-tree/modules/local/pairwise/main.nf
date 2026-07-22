process PAIRWISE {
    tag "${meta.id}"
    label 'process_high'

    conda "${projectDir}/environment.yml"

    input:
    tuple val(meta), path(chunk)

    output:
    tuple val(meta), path("${meta.id}.tmp.npz"), emit: matrices

    script:
    """
    Calculate_pairwise_dist_simple_large.py \
        --input ${chunk} \
        --output ${meta.id}.tmp.npz
    """

    stub:
    """
    touch ${meta.id}.tmp.npz
    """
}
