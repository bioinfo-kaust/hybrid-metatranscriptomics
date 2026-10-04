/*
 * Stage 1 — short-read QC: merge units (lanes) -> adapter/quality trim -> rRNA depletion -> FastQC.
 */

process CAT_FASTQ {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads1, stageAs: 'r1/*'), path(reads2, stageAs: 'r2/*')

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.merged.fq.gz"), emit: reads
    path "versions.yml", emit: versions, topic: versions

    script:
    // gzip members concatenate legally, so lane merging is a plain `cat`
    """
    cat r1/* > ${meta.id}_R1.merged.fq.gz
    cat r2/* > ${meta.id}_R2.merged.fq.gz
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch ${meta.id}_R1.merged.fq.gz ${meta.id}_R2.merged.fq.gz
    emit_versions.sh "${task.process}"
    """
}

process FASTP {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.trim.fq.gz"), emit: reads
    path "${meta.id}.fastp.json"                        , emit: json
    path "${meta.id}.fastp.html"                        , emit: html
    path "${meta.id}.fastp.log"                         , emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '--detect_adapter_for_pe'
    """
    fastp -i ${reads[0]} -I ${reads[1]} \\
          -o ${meta.id}_R1.trim.fq.gz -O ${meta.id}_R2.trim.fq.gz \\
          --thread ${task.cpus} ${args} \\
          --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html \\
          2> ${meta.id}.fastp.log
    emit_versions.sh "${task.process}" fastp
    """

    stub:
    """
    touch ${meta.id}_R1.trim.fq.gz ${meta.id}_R2.trim.fq.gz ${meta.id}.fastp.json ${meta.id}.fastp.html ${meta.id}.fastp.log
    emit_versions.sh "${task.process}" fastp
    """
}

/*
 * rRNA depletion with bbduk: a pair is dropped when EITHER mate matches an rRNA k-mer.
 * bbduk rather than SortMeRNA — native paired output, low memory, and it reports %rRNA per
 * database, which is what the QC report consumes.
 */
process BBDUK_RRNA_PE {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    path rrna_db

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.rrna_clean.fq.gz"), emit: reads
    path "${meta.id}.bbduk_stats.txt"                          , emit: stats
    path "${meta.id}.bbduk_refstats.txt"                       , emit: refstats
    path "${meta.id}.bbduk.log"                                , emit: log
    path "${meta.id}.rrna.tsv"                                 , emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    def mem  = (task.memory ? task.memory.toGiga() - 4 : 48)
    """
    bbduk.sh -Xmx${mem}g threads=${task.cpus} \\
        in1=${reads[0]} in2=${reads[1]} \\
        out1=${meta.id}_R1.rrna_clean.fq.gz out2=${meta.id}_R2.rrna_clean.fq.gz \\
        ref=${rrna_db} k=${params.rrna_kmer} ${args} \\
        stats=${meta.id}.bbduk_stats.txt refstats=${meta.id}.bbduk_refstats.txt \\
        overwrite=t 2> ${meta.id}.bbduk.log

    # per-sample rRNA percentage, parsed from bbduk's own accounting.
    # --paired: bbduk counts reads, every other per-sample table here counts pairs
    rrna_summary.py ${meta.id} ${meta.id}.bbduk.log ${meta.condition} --paired > ${meta.id}.rrna.tsv
    cat ${meta.id}.rrna.tsv
    emit_versions.sh "${task.process}" bbduk.sh
    """

    stub:
    """
    touch ${meta.id}_R1.rrna_clean.fq.gz ${meta.id}_R2.rrna_clean.fq.gz ${meta.id}.bbduk_stats.txt ${meta.id}.bbduk_refstats.txt ${meta.id}.bbduk.log
    printf "sample\tcondition\treads_in\trrna_removed\treads_kept\trrna_pct\n${meta.id}\t${meta.condition}\t1000\t100\t900\t10.00\n" > ${meta.id}.rrna.tsv
    emit_versions.sh "${task.process}" bbduk.sh
    """
}

/* Single-end / FASTA variant, used for the long-read isoforms. */
process BBDUK_RRNA_SE {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(seqs)
    path rrna_db

    output:
    tuple val(meta), path("${meta.id}.rrna_clean.fasta"), emit: seqs
    path "${meta.id}.rrna_hits.fasta"                   , emit: removed
    path "${meta.id}.bbduk_stats.txt"                   , emit: stats
    path "${meta.id}.bbduk.log"                         , emit: log
    path "${meta.id}.longread_rrna.tsv"                 , emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    def mem = (task.memory ? task.memory.toGiga() - 4 : 48)
    """
    bbduk.sh -Xmx${mem}g threads=${task.cpus} \\
        in=${seqs} out=${meta.id}.rrna_clean.fasta outm=${meta.id}.rrna_hits.fasta \\
        ref=${rrna_db} k=${params.rrna_kmer} \\
        stats=${meta.id}.bbduk_stats.txt refstats=${meta.id}.bbduk_refstats.txt \\
        overwrite=t 2> ${meta.id}.bbduk.log

    n_in=\$(grep -c '^>' ${seqs} || echo 0)
    n_rm=\$(grep -c '^>' ${meta.id}.rrna_hits.fasta || echo 0)
    n_ke=\$(grep -c '^>' ${meta.id}.rrna_clean.fasta || echo 0)
    printf 'sample\\tisoforms_in\\trrna_removed\\tkept\\tpct_removed\\n' > ${meta.id}.longread_rrna.tsv
    awk -v s=${meta.id} -v i=\$n_in -v r=\$n_rm -v k=\$n_ke \\
        'BEGIN{printf "%s\\t%d\\t%d\\t%d\\t%.2f\\n", s, i, r, k, (i? 100*r/i : 0)}' >> ${meta.id}.longread_rrna.tsv
    cat ${meta.id}.longread_rrna.tsv
    emit_versions.sh "${task.process}" bbduk.sh
    """

    stub:
    """
    touch ${meta.id}.rrna_clean.fasta ${meta.id}.rrna_hits.fasta ${meta.id}.bbduk_stats.txt ${meta.id}.bbduk.log
    printf "sample\tisoforms_in\trrna_removed\tkept\tpct_removed\n${meta.id}\t1000\t100\t900\t10.00\n" > ${meta.id}.longread_rrna.tsv
    emit_versions.sh "${task.process}" bbduk.sh
    """
}

process FASTQC {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)

    output:
    path '*_fastqc.{zip,html}', emit: report
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads}
    emit_versions.sh "${task.process}" fastqc
    """

    stub:
    """
    touch ${meta.id}_fastqc.zip ${meta.id}_fastqc.html
    emit_versions.sh "${task.process}" fastqc
    """
}

/*
 * Long reads arrive as FASTA or FASTQ (optionally gzipped) and may carry characters that break
 * downstream tools' shell handling (the project's own HQ file is named "...transcripts(1).fasta").
 * Normalise to a plain FASTA under a safe name before anything else touches it.
 */
process LONGREAD_PREP {
    label 'process_low'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fasta"), emit: seqs
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    src=\$(ls input/* | head -1)
    case "\$src" in
      *.gz)  zcat "\$src" > raw ;;
      *)     cp   "\$src"   raw ;;
    esac
    if head -1 raw | grep -q '^@'; then
        seqkit fq2fa raw > ${meta.id}.fasta
        rm -f raw
    else
        mv raw ${meta.id}.fasta
    fi
    grep -c '^>' ${meta.id}.fasta | sed 's/^/[${meta.id}] long reads: /'
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch ${meta.id}.fasta
    emit_versions.sh "${task.process}" seqkit
    """
}

/*
 * Which sequences came from which long-read sample. Long reads are pooled before they are mapped
 * (one pass over a multi-gigabase contaminant reference instead of one per sample), so per-sample
 * statistics are recovered afterwards from this map rather than from separate mapping runs.
 */
process LONGREAD_MEMBERSHIP {
    label 'process_single'
    tag "${meta.id}"

    input:
    tuple val(meta), path(seqs)

    output:
    path "${meta.id}.membership.tsv", emit: tsv
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    grep '^>' ${seqs} | sed 's/^>//; s/[[:space:]].*//' \\
      | awk -v s='${meta.id}' '{print \$1"\\t"s}' > ${meta.id}.membership.tsv
    echo "[${meta.id}] \$(wc -l < ${meta.id}.membership.tsv) sequences"
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    printf 'seq1\\t${meta.id}\\n' > ${meta.id}.membership.tsv
    emit_versions.sh "${task.process}" python3
    """
}

/* Pool the (usually few, unreplicated) long-read samples into one reference-building input. */
process LONGREAD_POOL {
    label 'process_low'

    input:
    path(fastas, stageAs: 'in/*')

    output:
    path 'longreads_pooled.fasta', emit: seqs
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    cat in/* > pooled.raw.fasta
    # duplicate names across pooled samples would silently collide downstream
    seqkit rmdup -n pooled.raw.fasta -o longreads_pooled.fasta 2> rmdup.log || mv pooled.raw.fasta longreads_pooled.fasta
    cat rmdup.log 2>/dev/null || true
    rm -f pooled.raw.fasta
    grep -c '^>' longreads_pooled.fasta | sed 's/^/[pooled] long reads: /'
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch longreads_pooled.fasta
    emit_versions.sh "${task.process}" seqkit
    """
}
