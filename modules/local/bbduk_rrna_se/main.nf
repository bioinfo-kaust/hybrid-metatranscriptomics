/*
 * rRNA depletion of long-read isoforms (FASTA, single-end). Kept local because the nf-core bbduk
 * module neither keeps the matched sequences (outm) nor takes FASTA output names; the removed set
 * feeds the per-sample long-read accounting.
 */
process BBDUK_RRNA_SE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bbmap_gawk:2749f656222129a6'
        : 'community.wave.seqera.io/library/bbmap_gawk:1f257e2d48ef74d8' }"

    input:
    tuple val(meta), path(seqs)
    path rrna_db

    output:
    tuple val(meta), path("${meta.id}.rrna_clean.fasta"), emit: seqs
    path "${meta.id}.rrna_hits.fasta"                   , emit: removed
    path "${meta.id}.bbduk_stats.txt"                   , emit: stats
    path "${meta.id}.bbduk.log"                         , emit: log
    path "${meta.id}.longread_rrna.tsv"                 , emit: summary
    path "versions.yml"                                 , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // leave the JVM headroom below the task's memory limit
    def mem  = task.memory ? Math.max(1, task.memory.toGiga() - 4) : 48
    """
    bbduk.sh -Xmx${mem}g threads=${task.cpus} \\
        in=${seqs} out=${meta.id}.rrna_clean.fasta outm=${meta.id}.rrna_hits.fasta \\
        ref=${rrna_db} ${args} \\
        stats=${meta.id}.bbduk_stats.txt refstats=${meta.id}.bbduk_refstats.txt \\
        overwrite=t 2> ${meta.id}.bbduk.log

    n_in=\$(grep -c '^>' ${seqs} || true)
    n_rm=\$(grep -c '^>' ${meta.id}.rrna_hits.fasta || true)
    n_ke=\$(grep -c '^>' ${meta.id}.rrna_clean.fasta || true)
    printf 'sample\\tisoforms_in\\trrna_removed\\tkept\\tpct_removed\\n' > ${meta.id}.longread_rrna.tsv
    awk -v s=${meta.id} -v i=\$n_in -v r=\$n_rm -v k=\$n_ke \\
        'BEGIN{printf "%s\\t%d\\t%d\\t%d\\t%.2f\\n", s, i, r, k, (i? 100*r/i : 0)}' >> ${meta.id}.longread_rrna.tsv
    cat ${meta.id}.longread_rrna.tsv
    emit_versions.sh "${task.process}" bbduk.sh
    """

    stub:
    """
    touch ${meta.id}.rrna_clean.fasta ${meta.id}.rrna_hits.fasta ${meta.id}.bbduk_stats.txt ${meta.id}.bbduk.log
    printf "sample\\tisoforms_in\\trrna_removed\\tkept\\tpct_removed\\n${meta.id}\\t1000\\t100\\t900\\t10.00\\n" > ${meta.id}.longread_rrna.tsv
    emit_versions.sh "${task.process}" bbduk.sh
    """
}
