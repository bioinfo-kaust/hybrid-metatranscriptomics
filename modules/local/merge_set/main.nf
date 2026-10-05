/* Concatenate a base set with the novel contributions of one or more sources. */
process MERGE_SET {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/seqkit_gawk:8a03f2b8f8f6d92a'
        : 'community.wave.seqera.io/library/seqkit_gawk:787aa427251d2eb0' }"

    input:
    tuple val(meta), path(parts, stageAs: 'part/*')

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: fasta
    path "versions.yml"                      , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // parts are concatenated in file-name order; vsearch centroids arrive gzipped
    """
    zcat -f part/* > ${meta.id}.fasta
    echo "[${meta.id}] transcripts: \$(grep -c '^>' ${meta.id}.fasta || true)"
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch ${meta.id}.fasta
    emit_versions.sh "${task.process}" seqkit
    """
}
