process FINALIZE_BQSR_BAM {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    tuple val(meta), path("${prefix}.bam.bai"), emit: bai

    script:
    prefix = task.ext.prefix ?: "${meta.id}_aln_sort_MD"
    """
    cp --reflink=auto ${bam} ${prefix}.bam
    cp --reflink=auto ${bai} ${prefix}.bam.bai
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_aln_sort_MD"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    """
}
