/*
 * Expression-filtered reference. The ladder's cross-source union carries every near-duplicate and
 * every unexpressed fragment into the index — which is what 31% BUSCO duplication looks like —
 * so the first quantification pass is used to drop what the experiment never saw, and the
 * survivors are re-indexed and re-quantified.
 */
process FILTER_BY_EXPRESSION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    path quants, stageAs: 'quants/*'
    tuple val(meta), path(fasta)
    tuple val(meta2), path(tx2gene)
    tuple val(meta3), path(annotation)

    output:
    tuple val(meta_out), path("${meta.id}_expressed.fasta")         , emit: fasta
    tuple val(meta_out), path("tx2gene.${meta.id}_expressed.tsv")   , emit: tx2gene
    tuple val(meta_out), path("annotation.${meta.id}_expressed.tsv"), emit: annotation
    path "${meta.id}_expressed.filter.txt"                          , emit: log
    path "${meta.id}_expressed.de_eligible_tx.txt"                  , emit: expressed_tx
    path "versions.yml"                                             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta_out  = meta + [id: "${meta.id}_expressed"]
    def label = meta.id
    def ann   = annotation.name != 'NO_FILE_ANN'
        ? "--annotation ${annotation} --out-annotation annotation.${label}_expressed.tsv" : ''
    // the rescue needs the annotation, so it silently degrades to expression-only without one
    def rescue = (params.expr_filter_rescue && annotation.name != 'NO_FILE_ANN')
        ? '--rescue-annotated' : ''
    """
    filter_by_expression.py \\
        --quant-dir quants --fasta ${fasta} --tx2gene ${tx2gene} ${ann} ${rescue} \\
        --min-tpm ${params.expr_filter_tpm} --min-samples ${params.expr_filter_samples} \\
        --out-fasta ${label}_expressed.fasta --out-tx2gene tx2gene.${label}_expressed.tsv \\
        --out-expressed-tx ${label}_expressed.de_eligible_tx.txt \\
        --log ${label}_expressed.filter.txt
    # a reference without annotation still needs the file to exist for the channel
    [ -s annotation.${label}_expressed.tsv ] || touch annotation.${label}_expressed.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    meta_out  = meta + [id: "${meta.id}_expressed"]
    def label = meta.id
    """
    touch ${label}_expressed.fasta tx2gene.${label}_expressed.tsv annotation.${label}_expressed.tsv ${label}_expressed.filter.txt ${label}_expressed.de_eligible_tx.txt
    emit_versions.sh "${task.process}" python3
    """
}
