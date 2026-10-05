/*
 * What this source ADDS to a reference set. Alignment-based rather than word-based: a query is
 * already represented when some alignment covers at least `novelty_min_cov` of it at the given
 * identity, and everything else is this source's contribution. minimap2 is already a dependency,
 * the rule is one readable threshold pair instead of cd-hit-est-2d's word-length heuristics, and
 * the per-query evidence stays in the PAF.
 */
process NOVEL_VS_REF {
    tag "${meta.id} vs ${meta2.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:49f77f39b14c7e4f'
        : 'community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:9b379af53c51bfba' }"

    input:
    tuple val(meta), path(query)
    tuple val(meta2), path(reference)
    val identity
    val prefix

    output:
    tuple val(meta), path("${meta.id}.novel.fasta"), emit: fasta
    path "${meta.id}.novelty.log"                  , emit: log
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: '-x asm20 -c --secondary=no'
    def label = meta.id
    // no reference yet (no genome-guided base): every query sequence is novel
    def align = reference.name == 'NO_FILE_REF' ? ': > q_vs_ref.paf'
              : "minimap2 ${args} -t ${task.cpus} ${reference} ${query} > q_vs_ref.paf 2> minimap2.log"
    """
    ${align}

    # PAF: 1 qname 2 qlen 3 qstart 4 qend .. 10 matching bases 11 alignment length
    awk -v cov=${params.novelty_min_cov} -v idt=${identity} \\
        '{ qcov = (\$4 - \$3) / \$2; qid = \$10 / \$11; if (qcov >= cov && qid >= idt) print \$1 }' \\
        q_vs_ref.paf | sort -u > represented.ids
    seqkit seq -n -i ${query} | sort -u > all.ids
    comm -23 all.ids represented.ids > novel.ids
    seqkit grep -f novel.ids ${query} > novel.raw.fasta

    # Sanitise the ID ONLY. Rewriting the whole header welds the assembler's description onto the
    # id (">TR_TRINITY_DN1_c0_g1_i1_len_328_path__0_0-327_"), and the tx2gene isoform rule is
    # anchored at the end of the id, so it then never matches and every isoform becomes its own
    # gene. The description carries no information the pipeline uses.
    awk -v p='${prefix}' '/^>/{id=substr(\$1,2); gsub(/[^A-Za-z0-9._-]/,"_",id); print ">" p id; next}{print}' novel.raw.fasta \\
      | seqkit seq -m ${params.min_transcript_len} - > ${label}.novel.fasta

    {
      echo "[${label}] vs ${meta2.id} at id>=${identity}, query coverage>=${params.novelty_min_cov}"
      echo "  query sequences      : \$(wc -l < all.ids)"
      echo "  already represented  : \$(wc -l < represented.ids)"
      echo "  novel (>=${params.min_transcript_len} bp): \$(grep -c '^>' ${label}.novel.fasta || true)"
    } | tee ${label}.novelty.log
    emit_versions.sh "${task.process}" minimap2 seqkit
    """

    stub:
    """
    touch ${meta.id}.novel.fasta ${meta.id}.novelty.log
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}
