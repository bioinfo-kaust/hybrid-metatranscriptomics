/* A directory of contaminant FASTA files (plain or gzipped), combined as-is. */
process COMBINE_CONTAMINANTS {
    tag "${fasta_dir}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/seqkit_gawk:8a03f2b8f8f6d92a'
        : 'community.wave.seqera.io/library/seqkit_gawk:787aa427251d2eb0' }"

    input:
    path fasta_dir

    output:
    path 'contaminants_combined.fna', emit: fasta
    path "versions.yml"             , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    shopt -s nullglob
    files=( ${fasta_dir}/*.fa ${fasta_dir}/*.fna ${fasta_dir}/*.fasta ${fasta_dir}/*.fa.gz ${fasta_dir}/*.fna.gz ${fasta_dir}/*.fasta.gz )
    [ \${#files[@]} -gt 0 ] || { echo "[ERROR] no FASTA files in ${fasta_dir}" >&2; exit 1; }
    : > contaminants_combined.fna
    for f in "\${files[@]}"; do
        zcat -f "\$f" >> contaminants_combined.fna
    done
    echo "[contaminants] scaffolds: \$(grep -c '^>' contaminants_combined.fna || true)"
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch contaminants_combined.fna
    emit_versions.sh "${task.process}" seqkit
    """
}
