process CREATE_SUBSET_SAMPLE_ARGS {
    tag 'subset'
    label 'process_single'

    input:
    val sample_ids

    output:
    path 'subset_samples.args', emit: sample_args

    script:
    def rows = sample_ids.join('\n') + '\n'
    def encoded_rows = rows.bytes.encodeBase64().toString()
    """
    printf '%s' '${encoded_rows}' | base64 --decode > subset_samples.args
    """

    stub:
    """
    touch subset_samples.args
    """
}
