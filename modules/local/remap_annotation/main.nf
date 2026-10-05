/*
 * Re-key an annotation table onto a different transcript->gene map.
 *
 * Its own process on purpose: folding this into CORSET made corset's 37-minute clustering
 * re-run every time the annotation or the remap changed, for a result corset does not even read.
 * Anything downstream of a gene map can now change without touching the expensive step.
 */
process REMAP_ANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(tx2gene)
    tuple val(meta2), path(annotation)

    output:
    tuple val(meta), path("annotation.${meta.id}.tsv"), emit: annotation
    path "versions.yml"                               , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label = meta.id
    """
    # The enrichment keys on the annotation's gene column while the DE table is keyed on the new
    # map's ids: leave them unmatched and the universe comes out empty.
    awk -F'\\t' -v OFS='\\t' '
        NR==FNR { g[\$1] = \$2; next }
        FNR==1  { print; next }
        { if (\$1 in g) \$2 = g[\$1]; print }
    ' ${tx2gene} ${annotation} > annotation.${label}.tsv

    echo "[remap] ${label}: \$(( \$(wc -l < annotation.${label}.tsv) - 1 )) rows re-keyed onto \$(cut -f2 ${tx2gene} | sort -u | wc -l) genes"
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch annotation.${meta.id}.tsv
    emit_versions.sh "${task.process}"
    """
}
