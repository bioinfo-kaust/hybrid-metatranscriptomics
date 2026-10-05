/*
 * Route A models with GENOME bases: gffread spliced at the collapsed models' exon coordinates.
 * The Iso-Seq consensus carries homopolymer indels that shift the reading frame; the genome at
 * the same coordinates does not (this is what SQANTI3 writes as *_corrected.fasta). Ids and model
 * structure are unchanged, so every downstream map stays valid. Fails if any model is lost.
 */
process LONGREAD_GENOME_SEQ {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/gffread_seqkit_gawk:97e9b64ec67d60ea'
        : 'community.wave.seqera.io/library/gffread_seqkit_gawk:7d0bb4806e001c3b' }"

    input:
    tuple val(meta), path(gtf)
    path genome

    output:
    tuple val(meta), path('collapsed.genome.fasta'), emit: fasta
    path 'collapsed.genome.stats.tsv'              , emit: stats
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gffread -w raw.fa -g ${genome} ${gtf}
    # bare ids (gffread may append attributes) and upper case (soft-masked genomes give lower case)
    seqkit seq -i -u raw.fa > collapsed.genome.fasta
    n_models=\$(awk '\$3=="transcript"' ${gtf} | wc -l)
    n_seq=\$(grep -c '^>' collapsed.genome.fasta || true)
    [ "\$n_models" = "\$n_seq" ] || { echo "[ERROR] \$n_models models, \$n_seq genome sequences" >&2; exit 1; }
    seqkit stats -T collapsed.genome.fasta > collapsed.genome.stats.tsv
    cat collapsed.genome.stats.tsv
    emit_versions.sh "${task.process}" gffread seqkit
    """

    stub:
    """
    touch collapsed.genome.fasta collapsed.genome.stats.tsv
    emit_versions.sh "${task.process}" gffread seqkit
    """
}
