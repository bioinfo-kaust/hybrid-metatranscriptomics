/* A protein FASTA for DIAMOND (Swiss-Prot, a specialist set, or anything else). */
process DOWNLOAD_PROTEIN_DB {
    tag "${meta.id}"
    label 'process_single'
    label 'process_download'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/curl_gzip:a0a9da312d0da956'
        : 'community.wave.seqera.io/library/curl_gzip:7d4ca65cc6f11ae1' }"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("${meta.id}.fasta.gz"), emit: fasta
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    curl -sfL --max-time 1800 "${url}" -o ${meta.id}.download
    # UniProt streams can arrive uncompressed despite compressed=true
    if gzip -t ${meta.id}.download 2>/dev/null; then
        mv ${meta.id}.download ${meta.id}.fasta.gz
    else
        gzip -c ${meta.id}.download > ${meta.id}.fasta.gz && rm ${meta.id}.download
    fi
    [ -s ${meta.id}.fasta.gz ] || { echo "[ERROR] empty download for ${meta.id}" >&2; exit 1; }
    emit_versions.sh "${task.process}" curl
    """

    stub:
    """
    echo "" | gzip > ${meta.id}.fasta.gz
    emit_versions.sh "${task.process}" curl
    """
}
