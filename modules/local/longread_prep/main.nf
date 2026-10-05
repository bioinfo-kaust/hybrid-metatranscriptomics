/*
 * Long reads arrive as FASTA or FASTQ (optionally gzipped) and may carry characters that break
 * downstream tools' shell handling (the project's own HQ file is named "...transcripts(1).fasta").
 * Normalise to a plain FASTA under a safe name before anything else touches it.
 */
process LONGREAD_PREP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/seqkit_gawk:8a03f2b8f8f6d92a'
        : 'community.wave.seqera.io/library/seqkit_gawk:787aa427251d2eb0' }"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: seqs
    path "versions.yml"                      , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    files=( input/* ); src="\${files[0]}"
    zcat -f "\$src" > raw
    if head -c 1 raw | grep -q '@'; then
        seqkit fq2fa raw > ${meta.id}.fasta
        rm -f raw
    else
        mv raw ${meta.id}.fasta
    fi
    echo "[${meta.id}] long reads: \$(grep -c '^>' ${meta.id}.fasta || true)"
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch ${meta.id}.fasta
    emit_versions.sh "${task.process}" seqkit
    """
}
