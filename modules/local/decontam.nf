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

    input:
    path genome
    path contaminants

    output:
    path 'host_plus_contam.fna'   , emit: fasta
    path 'host_plus_contam.sr.mmi', emit: index
    path 'decon_index.log'        , emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    def copy_host = genome.name == 'NO_FILE_GENOME' ? ': > host_plus_contam.fna' : "cp ${genome} host_plus_contam.fna"
    """
    ${copy_host}
    awk '/^>/{sub(/^>/,">${params.contaminant_prefix}"); print; next}{print}' ${contaminants} >> host_plus_contam.fna
    {
      echo "scaffolds total:       \$(grep -c '^>' host_plus_contam.fna)"
      echo "scaffolds contaminant: \$(grep -c '^>${params.contaminant_prefix}' host_plus_contam.fna)"
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

process DECONTAM_SHORT {
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.host.fq.gz"), emit: reads
    path "${meta.id}.decontam.tsv"                       , emit: summary
    path "${meta.id}.minimap2.log"                       , emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    # -F 0x900 keeps primary, non-supplementary records: one best placement per read
    minimap2 -ax sr -t ${task.cpus} ${index} ${reads[0]} ${reads[1]} 2> ${meta.id}.minimap2.log \\
      | samtools view -F 0x900 - \\
      | awk -F'\\t' '\$3 ~ /^${params.contaminant_prefix}/ {print \$1}' \\
      | LC_ALL=C sort -u -S 4G -T . > contaminant_read_names.txt
    n_contam=\$(wc -l < contaminant_read_names.txt)
    echo "[${meta.id}] read pairs placed on a contaminant: \$n_contam"

    for m in 1 2; do
        src=\$([ \$m = 1 ] && echo ${reads[0]} || echo ${reads[1]})
        if [ "\$n_contam" -gt 0 ]; then
            seqkit grep -v -f contaminant_read_names.txt "\$src" -o ${meta.id}_R\${m}.host.fq.gz
        else
            cp "\$src" ${meta.id}_R\${m}.host.fq.gz
        fi
    done

    nseq () { seqkit stats -T "\$1" | awk 'NR==2{print \$4}'; }
    n_in=\$(nseq ${reads[0]})
    n_k1=\$(nseq ${meta.id}_R1.host.fq.gz)
    n_k2=\$(nseq ${meta.id}_R2.host.fq.gz)
    [ "\$n_k1" = "\$n_k2" ] || { echo "[ERROR] ${meta.id}: kept R1/R2 counts differ (\$n_k1 vs \$n_k2) — pairing broken" >&2; exit 1; }

    printf 'sample\\tcondition\\tpairs_in\\tcontaminant_removed\\tpairs_kept\\tpct_removed\\n' > ${meta.id}.decontam.tsv
    awk -v s=${meta.id} -v c=${meta.condition} -v i=\$n_in -v r=\$n_contam -v k=\$n_k1 \\
        'BEGIN{printf "%s\\t%s\\t%d\\t%d\\t%d\\t%.2f\\n", s, c, i, r, k, (i? 100*r/i : 0)}' >> ${meta.id}.decontam.tsv
    cat ${meta.id}.decontam.tsv
    emit_versions.sh "${task.process}" minimap2 seqkit
    """

    stub:
    """
    touch ${meta.id}_R1.host.fq.gz ${meta.id}_R2.host.fq.gz ${meta.id}.minimap2.log
    printf "sample\tcondition\tpairs_in\tcontaminant_removed\tpairs_kept\tpct_removed\n${meta.id}\t${meta.condition}\t1000\t130\t870\t13.00\n" > ${meta.id}.decontam.tsv
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}

/*
 * Long reads: spliced alignment to the host genome and to the contaminant set, then the
 * best-residue-match classifier. 'contaminant' (and, by default, 'ambiguous') are dropped;
 * 'host' and unaligned are kept — taking the complement matters, because the classifier only
 * enumerates sequences that appear in at least one PAF.
 */
process DECONTAM_LONG {
    label 'process_high'
    tag "${label}"

    input:
    tuple val(label), path(seqs)
    path genome
    path contaminants
    path membership, stageAs: 'membership/*'
    path rrna_hits , stageAs: 'rrna_hits/*'
    path conditions

    output:
    tuple val(label), path("${label}.hostonly.fasta"), emit: seqs
    path "${label}.contaminant_removed.fasta"        , emit: removed
    path "${label}.longread_decontam.tsv"            , emit: summary
    path "${label}.classification.tsv"               , emit: classification
    path "longread_decontam_per_sample.tsv"          , emit: per_sample, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    def drop_ambiguous = params.keep_ambiguous_long ? '' : "${label}.ambiguous.ids"
    def keep_amb  = params.keep_ambiguous_long ? '--keep-ambiguous' : ''
    def cond_arg  = conditions.name != 'NO_FILE_COND' ? "--conditions ${conditions}" : ''
    // no genome: an empty host PAF, so every contaminant placement classifies as 'contaminant'
    def host_aln  = genome.name == 'NO_FILE_GENOME' ? ': > vs_host.paf'
                  : "minimap2 -x splice -t ${task.cpus} --secondary=no ${genome} ${seqs} > vs_host.paf 2> vs_host.log"
    """
    ${host_aln}
    minimap2 -x splice -t ${task.cpus} --secondary=no ${contaminants} ${seqs} > vs_contam.paf 2> vs_contam.log

    classify_host_contaminant.py vs_host.paf vs_contam.paf ${label} ${params.contaminant_margin}

    cat ${label}.contaminant.ids ${drop_ambiguous} 2>/dev/null | sort -u > ${label}.removed.ids
    if [ -s ${label}.removed.ids ]; then
        seqkit grep -v -f ${label}.removed.ids ${seqs} > ${label}.hostonly.fasta
        seqkit grep    -f ${label}.removed.ids ${seqs} > ${label}.contaminant_removed.fasta
    else
        cp ${seqs} ${label}.hostonly.fasta
        : > ${label}.contaminant_removed.fasta
    fi

    N_IN=\$(grep -c '^>' ${seqs})
    N_H=\$(wc -l < ${label}.host.ids); N_S=\$(wc -l < ${label}.contaminant.ids); N_A=\$(wc -l < ${label}.ambiguous.ids)
    {
      printf 'class\\tn_sequences\\n'
      printf 'host\\t%s\\n'            "\$N_H"
      printf 'no_alignment\\t%s\\n'    "\$(( N_IN - N_H - N_S - N_A ))"
      printf 'contaminant\\t%s\\n'     "\$N_S"
      printf 'ambiguous\\t%s\\n'       "\$N_A"
      printf 'kept\\t%s\\n'            "\$(grep -c '^>' ${label}.hostonly.fasta)"
      printf 'removed\\t%s\\n'         "\$(grep -c '^>' ${label}.contaminant_removed.fasta || echo 0)"
    } | tee ${label}.longread_decontam.tsv

    # per-sample accounting, recovered from the membership map rather than from separate runs
    if ls membership/* >/dev/null 2>&1; then
        longread_per_sample.py --membership membership --classification ${label}.classification.tsv \\
            --input-ids ${seqs} --rrna-dir rrna_hits ${cond_arg} ${keep_amb} \\
            --out longread_decontam_per_sample.tsv
    fi
    emit_versions.sh "${task.process}" minimap2 seqkit
    """

    stub:
    """
    touch ${label}.hostonly.fasta ${label}.contaminant_removed.fasta ${label}.classification.tsv
    printf "class\tn_sequences\nhost\t100\n" > ${label}.longread_decontam.tsv
    printf "sample\tcondition\tisoforms_in\trrna_removed\thost\tno_alignment\tcontaminant\tambiguous\tcontaminant_removed\tkept\tpct_removed\n" > longread_decontam_per_sample.tsv
    printf "LR1\tcond\t1000\t10\t800\t100\t80\t10\t90\t900\t9.09\n" >> longread_decontam_per_sample.tsv
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}

/*
 * Residual audit of an assembled set. After read-level removal this is a QC gate, not the
 * separation step: it reports how much contaminant signal survived, and drops what did.
 * keep_class controls the keep rule ('host,ambiguous' matches the project's audit).
 */
process CLASSIFY_CONTAM {
    label 'process_high'
    tag "$label"

    input:
    tuple val(label), path(fasta)
    path genome
    path contaminants
    val keep_class

    output:
    tuple val(label), path("${label}.host.fasta"), emit: fasta
    path "${label}.hostsym_counts.txt"           , emit: counts
    path "${label}.classification.tsv"           , emit: classification
    path "${label}.host.ids"                     , emit: host_ids
    path "versions.yml", emit: versions, topic: versions

    script:
    def host_aln = genome.name == 'NO_FILE_GENOME' ? ': > vs_host.paf'
                 : "minimap2 -x splice -t ${task.cpus} --secondary=no ${genome} ${fasta} > vs_host.paf 2> vs_host.log"
    """
    ${host_aln}
    minimap2 -x splice -t ${task.cpus} --secondary=no ${contaminants} ${fasta} > vs_contam.paf 2> vs_contam.log
    classify_host_contaminant.py vs_host.paf vs_contam.paf ${label} ${params.contaminant_margin}

    grep '^>' ${fasta} | sed 's/^>//; s/[[:space:]].*//' | sort -u > all.ids

    # A sequence that aligns to NEITHER reference never appears in either PAF, so the classifier
    # never lists it in any class — it is not "unmapped: 0", it is absent. Derive that class from
    # the FASTA instead, or the keep rule deletes it: 19,166 of 129,450 Trinity contigs (14.8%)
    # and 2,804 of 79,140 reference-free long-read sequences on this dataset, which is exactly the
    # genome-absent material the de novo route exists to recover.
    cat ${label}.host.ids ${label}.contaminant.ids ${label}.ambiguous.ids 2>/dev/null \\
      | sort -u > classified.ids
    comm -23 all.ids classified.ids > ${label}.unmapped.ids

    : > keep.ids
    for c in \$(echo "${keep_class}" | tr ',' ' '); do
        # `[ -s f ] && cat f` as the last command in the body exits non-zero on an empty class
        # file, which under `bash -ue` aborts the task
        if [ -f "${label}.\$c.ids" ]; then cat "${label}.\$c.ids" >> keep.ids; fi
    done
    sort -u keep.ids -o keep.ids
    seqkit grep -f keep.ids ${fasta} > ${label}.host.fasta

    {
      echo "[${label}] classification (post read-level removal)"
      for c in host contaminant ambiguous unmapped; do
          echo "  \$c: \$(wc -l < ${label}.\$c.ids)"
      done
      echo "  input sequences: \$(wc -l < all.ids)"
      echo "  kept (${keep_class}): \$(grep -c '^>' ${label}.host.fasta)"
      echo "  dropped: \$(( \$(wc -l < all.ids) - \$(grep -c '^>' ${label}.host.fasta) ))"
    } | tee ${label}.hostsym_counts.txt
    emit_versions.sh "${task.process}" minimap2 seqkit
    """

    stub:
    """
    touch ${label}.host.fasta ${label}.hostsym_counts.txt ${label}.classification.tsv ${label}.host.ids
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}
