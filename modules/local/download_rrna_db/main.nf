/*
 * rRNA reference for k-mer depletion. Default set = SortMeRNA v4.3 (SILVA 16S/18S/23S/28S +
 * Rfam 5S/5.8S); any other set of URLs works via --rrna_db_sets / --rrna_db_base.
 */
process DOWNLOAD_RRNA_DB {
    label 'process_single'
    label 'process_download'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/curl_gzip:a0a9da312d0da956'
        : 'community.wave.seqera.io/library/curl_gzip:7d4ca65cc6f11ae1' }"

    input:
    val base_url
    val sets

    output:
    path 'rrna_db.fasta', emit: fasta
    path "versions.yml" , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    : > rrna_db.fasta
    for f in \$(echo "${sets}" | tr ',' ' '); do
        echo "[*] \$f"
        curl -sfL --max-time 600 "${base_url}/\$f.fasta" >> rrna_db.fasta
    done
    [ -s rrna_db.fasta ] || { echo "[ERROR] rRNA database download produced an empty file" >&2; exit 1; }
    echo "[rRNA db] sequences: \$(grep -c '^>' rrna_db.fasta || true)"
    emit_versions.sh "${task.process}" curl
    """

    stub:
    """
    touch rrna_db.fasta
    emit_versions.sh "${task.process}" curl
    """
}
