/*
 * Contaminant genomes to subtract, from a comma-separated accession list. The genomes are
 * downloaded and concatenated; BUILD_DECON_INDEX prefixes their scaffolds so a contaminant hit is
 * a prefix test.
 */
process DOWNLOAD_CONTAMINANT_GENOMES {
    tag "${accessions}"
    label 'process_single'
    label 'process_download'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/ncbi-datasets-cli_unzip_gawk:b4a8163c7e5795ed'
        : 'community.wave.seqera.io/library/ncbi-datasets-cli_unzip_gawk:bd7c22975b1f1721' }"

    input:
    val accessions

    output:
    path 'contaminants_combined.fna', emit: fasta
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    : > contaminants_combined.fna
    for acc in \$(echo "${accessions}" | tr ',' ' '); do
        echo "[*] datasets download \$acc"
        datasets download genome accession "\$acc" --include genome --filename "\$acc.zip"
        unzip -o -q "\$acc.zip" -d "\$acc"_dir
        cat "\$acc"_dir/ncbi_dataset/data/"\$acc"/*_genomic.fna >> contaminants_combined.fna
        rm -rf "\$acc.zip" "\$acc"_dir
    done
    echo "[contaminants] scaffolds: \$(grep -c '^>' contaminants_combined.fna || true)"
    emit_versions.sh "${task.process}" datasets
    """

    stub:
    """
    touch contaminants_combined.fna
    emit_versions.sh "${task.process}" datasets
    """
}
