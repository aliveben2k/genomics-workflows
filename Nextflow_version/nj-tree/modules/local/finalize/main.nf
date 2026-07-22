process FINALIZE {
    tag "${output_prefix}"
    label 'process_high'

    conda "${projectDir}/environment.yml"

    input:
    path matrices
    val output_prefix

    output:
    path "${output_prefix}.mtx", emit: distance_matrix
    path "${output_prefix}_final.npz", emit: archive
    path "${output_prefix}_pcoa.txt", emit: pcoa_table
    path "${output_prefix}_pcoa.pdf", emit: pcoa_plot
    path "${output_prefix}_tree.newick", emit: newick
    path "${output_prefix}_tree.nex", emit: nexus
    path "${output_prefix}_tree_pcoa.pkl", emit: python_results

    script:
    def matrix_args = matrices.collect { it.toString() }.join(' ')
    """
    vcf2table_large.py finalize-files \
        --output ${output_prefix} \
        --matrices ${matrix_args}
    """

    stub:
    """
    touch ${output_prefix}.mtx
    touch ${output_prefix}_final.npz
    touch ${output_prefix}_pcoa.txt
    touch ${output_prefix}_pcoa.pdf
    touch ${output_prefix}_tree.newick
    touch ${output_prefix}_tree.nex
    touch ${output_prefix}_tree_pcoa.pkl
    """
}
