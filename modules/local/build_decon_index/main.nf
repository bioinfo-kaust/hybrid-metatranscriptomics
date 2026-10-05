/*
 * Stage 2 — remove contaminant (e.g. endosymbiont) reads BEFORE anything is assembled.
 *
 * The reference genome and the contaminant genomes go into ONE index, so each read's primary
 * alignment is its best placement across both. A conserved host read that also aligns to a
 * contaminant still scores best on the host and is kept; only a positive contaminant placement
 * removes a pair. Reads that align nowhere are kept — on a fragmented draft genome,
 * unmapped is not evidence of contamination.
 *
 * Without a genome (genome input = the NO_FILE_GENOME placeholder) the index holds only the
 * contaminants and the host alignments are empty, so ANY placement on a contaminant removes the
 * read. There is no host to win a conserved region back; that is the price of running genome-free.
 */
process BUILD_DECON_INDEX {
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:49f77f39b14c7e4f'
        : 'community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:9b379af53c51bfba' }"

    input:
    path genome
    path contaminants
    val contaminant_prefix

    output:
    path 'host_plus_contam.fna'   , emit: fasta
    path 'host_plus_contam.sr.mmi', emit: index
    path 'decon_index.log'        , emit: log
    path "versions.yml"           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def copy_host = genome.name == 'NO_FILE_GENOME' ? ': > host_plus_contam.fna' : "zcat -f ${genome} > host_plus_contam.fna"
    """
    ${copy_host}
    zcat -f ${contaminants} | awk '/^>/{sub(/^>/,">${contaminant_prefix}"); print; next}{print}' >> host_plus_contam.fna
    {
      echo "scaffolds total:       \$(grep -c '^>' host_plus_contam.fna || true)"
      echo "scaffolds contaminant: \$(grep -c '^>${contaminant_prefix}' host_plus_contam.fna || true)"
    } | tee decon_index.log
    minimap2 -x sr -t ${task.cpus} -d host_plus_contam.sr.mmi host_plus_contam.fna
    emit_versions.sh "${task.process}" minimap2
    """

    stub:
    """
    touch host_plus_contam.fna host_plus_contam.sr.mmi decon_index.log
    emit_versions.sh "${task.process}" minimap2
    """
}
