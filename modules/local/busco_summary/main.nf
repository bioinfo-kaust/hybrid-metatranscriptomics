/*
 * The one-line completeness string from a BUSCO short summary, as <label>.busco.txt — the form the
 * ladder table, the figures and the dashboard read.
 */
process BUSCO_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(summary, stageAs: 'summary/*')

    output:
    tuple val(meta), path("${meta.id}.busco.txt"), emit: summary
    path "versions.yml"                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    { grep -hoE "C:[0-9.]+%\\[S:[0-9.]+%,D:[0-9.]+%\\],F:[0-9.]+%,M:[0-9.]+%" summary/*.txt 2>/dev/null || true; } \\
      | sed -n 1p > ${meta.id}.busco.txt
    [ -s ${meta.id}.busco.txt ] || echo NA >> ${meta.id}.busco.txt
    echo "[${meta.id}] \$(cat ${meta.id}.busco.txt)"
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    echo "C:70.0%[S:50.0%,D:20.0%],F:10.0%,M:20.0%" > ${meta.id}.busco.txt
    emit_versions.sh "${task.process}"
    """
}
