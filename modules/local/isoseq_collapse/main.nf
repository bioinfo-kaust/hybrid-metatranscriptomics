/*
 * Route A: collapse the genome-aligned long reads into isoform models.
 */
process ISOSEQ_COLLAPSE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/isoseq_seqkit_gawk:d6d0594b2bc0361c'
        : 'community.wave.seqera.io/library/isoseq_seqkit_gawk:9c63db4dfee47a7d' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('collapsed.gff')      , emit: gff
    tuple val(meta), path('collapsed.gtf')      , emit: gtf
    tuple val(meta), path('collapsed.fasta')    , emit: fasta
    path 'collapsed.stats.tsv'                  , emit: stats
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    isoseq collapse ${args} ${bam} collapsed.gff
    # IsoformSwitchAnalyzeR wants a .gtf; isoseq already writes GFF2-style attributes
    cp collapsed.gff collapsed.gtf
    # `isoseq collapse` decorates headers (PB.1.1|PB.1.1:5118-6219(-)|transcript/569365).
    # Reduce to the bare accession: the gene id is derived by stripping the trailing isoform
    # field, which only works on the bare form — leaving the decoration makes every isoform
    # its own gene and silently turns gene-level DGE into isoform-level.
    awk '/^>/{split(substr(\$0,2),a,"|"); print ">"a[1]; next}{print}' collapsed.fasta > tmp.fa
    mv tmp.fa collapsed.fasta
    seqkit stats -T collapsed.fasta > collapsed.stats.tsv
    cat collapsed.stats.tsv
    emit_versions.sh "${task.process}" isoseq seqkit
    """

    stub:
    """
    touch collapsed.gff collapsed.gtf collapsed.fasta collapsed.stats.tsv
    emit_versions.sh "${task.process}" isoseq seqkit
    """
}
