/*
 * The analysed reference, collapsed to one sequence per gene.
 *
 * A gene map is a table of ids: after it, the gene exists as a sum of counts and never as a
 * sequence. This writes the gene-level FASTA beside the map that defines it, so a gene set can be
 * BLASTed, scored by BUSCO in protein mode or shipped as a supplement without anyone rebuilding
 * it from the transcript FASTA by hand. The representative is the longest transcript (and the
 * longest ORF on any of the gene's transcripts for the peptide) — see gene_representatives.py for
 * why, and for what a representative is not.
 *
 * Run once per gene map, so the file always carries the map's own label in its name.
 */
process GENE_FASTA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(tx2gene)
    tuple val(meta2), path(fasta)
    path pep

    output:
    tuple val(meta), path("${meta.id}.genes.fasta")   , emit: fasta
    tuple val(meta), path("${meta.id}.genes.pep")     , emit: pep
    tuple val(meta), path("${meta.id}.gene_rep.tsv")  , emit: table
    path "versions.yml"                               , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label = meta.id
    """
    gene_representatives.py \\
        --fasta ${fasta} \\
        --tx2gene ${tx2gene} \\
        --pep ${pep} \\
        --out-fasta ${label}.genes.fasta \\
        --out-pep ${label}.genes.pep \\
        --out-table ${label}.gene_rep.tsv
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    printf ">LOCUS_000001 rep_transcript=t1 len=4 n_tx=1\\nACGT\\n" > ${meta.id}.genes.fasta
    touch ${meta.id}.genes.pep
    printf "gene\\tn_transcripts\\trep_transcript\\trep_tx_len\\trep_orf\\trep_orf_aa\\n" > ${meta.id}.gene_rep.tsv
    emit_versions.sh "${task.process}"
    """
}
