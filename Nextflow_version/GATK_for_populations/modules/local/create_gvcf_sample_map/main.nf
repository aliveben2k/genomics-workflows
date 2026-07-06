process CREATE_GVCF_SAMPLE_MAP {
    tag 'cohort'
    label 'process_single'

    input:
    tuple val(sample_ids), path(gvcfs)

    output:
    path 'cohort.sample_map', emit: sample_map

    script:
    def rows = sample_ids.withIndex().collect { sample_id, index ->
        "${sample_id}\t${gvcfs[index].name}"
    }.join('\n')
    def encoded_rows = rows.bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded_rows}' | base64 --decode > cohort.sample_map
    """

    stub:
    """
    touch cohort.sample_map
    """
}
