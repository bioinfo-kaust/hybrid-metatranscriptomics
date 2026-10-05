/*
 * Which sequences came from which long-read sample. Long reads are pooled before they are mapped
 * (one pass over a multi-gigabase contaminant reference instead of one per sample), so per-sample
 * statistics are recovered afterwards from this map rather than from separate mapping runs.
 */
process LONGREAD_MEMBERSHIP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(seqs)

    output:
    path "${meta.id}.membership.tsv", emit: tsv
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    { grep '^>' ${seqs} || true; } | sed 's/^>//; s/[[:space:]].*//' \\
      | awk -v s='${meta.id}' '{print \$1"\\t"s}' > ${meta.id}.membership.tsv
    echo "[${meta.id}] \$(wc -l < ${meta.id}.membership.tsv) sequences"
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    printf 'seq1\\t${meta.id}\\n' > ${meta.id}.membership.tsv
    emit_versions.sh "${task.process}" python3
    """
}
