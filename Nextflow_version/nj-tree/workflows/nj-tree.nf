include { PREPARE  } from '../modules/local/prepare/main'
include { PAIRWISE } from '../modules/local/pairwise/main'
include { FINALIZE } from '../modules/local/finalize/main'


workflow NJ_TREE {
    take:
    ch_input
    output_prefix
    window_size

    main:
    PREPARE(ch_input, output_prefix, window_size)

    ch_chunks = PREPARE.out.chunks.flatMap { meta, chunks ->
        def chunk_list = chunks instanceof List ? chunks : [chunks]
        chunk_list.collect { chunk ->
            def chunk_id = chunk.name.replaceFirst(/\.tmp\.txt\.gz$/, '')
            tuple(meta + [id: chunk_id], chunk)
        }
    }

    PAIRWISE(ch_chunks)

    ch_matrices = PAIRWISE.out.matrices
        .map { meta, matrix -> matrix }
        .collect()

    FINALIZE(ch_matrices, output_prefix)

    emit:
    chunks          = PREPARE.out.chunks
    manifest        = PREPARE.out.manifest
    chunk_matrices  = PAIRWISE.out.matrices
    distance_matrix = FINALIZE.out.distance_matrix
    archive         = FINALIZE.out.archive
    pcoa_table      = FINALIZE.out.pcoa_table
    pcoa_plot       = FINALIZE.out.pcoa_plot
    newick          = FINALIZE.out.newick
    nexus           = FINALIZE.out.nexus
    python_results  = FINALIZE.out.python_results
}
