/*
 * Host genome from NCBI by accession. Nothing here is species-specific: every input is a
 * parameter, and anything that can be fetched from a public source is fetched rather than
 * assumed to be staged already.
 */
process DOWNLOAD_GENOME {
    tag "${accession}"
    label 'process_single'
    label 'process_download'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/ncbi-datasets-cli_unzip_gawk:b4a8163c7e5795ed'
        : 'community.wave.seqera.io/library/ncbi-datasets-cli_unzip_gawk:bd7c22975b1f1721' }"

    input:
    val accession
    val name

    output:
    path "${name}.fna"  , emit: fasta
    path "versions.yml" , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    datasets download genome accession ${accession} --include genome --filename ${accession}.zip
    unzip -o -q ${accession}.zip -d ${accession}_dir
    cat ${accession}_dir/ncbi_dataset/data/${accession}/*_genomic.fna > ${name}.fna
    rm -rf ${accession}.zip ${accession}_dir
    echo "[${name}] scaffolds: \$(grep -c '^>' ${name}.fna || true)"
    emit_versions.sh "${task.process}" datasets
    """

    stub:
    """
    touch ${name}.fna
    emit_versions.sh "${task.process}" datasets
    """
}
