process GATK4_SELECTVARIANTS_FROM_GENOMICSDB {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(genomicsdb)
    path sample_args
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.g.vcf.gz"), emit: gvcf
    tuple val(meta), path("*.g.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.subset"
    def avail_mem = task.memory ? (task.memory.mega * 0.8).intValue() : 3072
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        SelectVariants \\
        --variant gendb://${genomicsdb} \\
        --sample-name ${sample_args} \\
        --reference ${fasta} \\
        --output ${prefix}.g.vcf.gz \\
        --create-output-variant-index true \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.subset"
    """
    echo | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi
    """
}
