/*
 * Long reads: spliced alignment to the host genome and to the contaminant set, then the
 * best-residue-match classifier. 'contaminant' (and, by default, 'ambiguous') are dropped;
 * 'host' and unaligned are kept — taking the complement matters, because the classifier only
 * enumerates sequences that appear in at least one PAF.
 */
process DECONTAM_LONG {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:49f77f39b14c7e4f'
        : 'community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:9b379af53c51bfba' }"

    input:
    tuple val(meta), path(seqs)
    path genome
    path contaminants
    path membership, stageAs: 'membership/*'
    path rrna_hits , stageAs: 'rrna_hits/*'
    path conditions

    output:
    tuple val(meta), path("${meta.id}.hostonly.fasta"), emit: seqs
    path "${meta.id}.contaminant_removed.fasta"       , emit: removed
    path "${meta.id}.longread_decontam.tsv"           , emit: summary
    path "${meta.id}.classification.tsv"              , emit: classification
    path "longread_decontam_per_sample.tsv"           , emit: per_sample, optional: true
    path "versions.yml"                               , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label          = meta.id
    def drop_ambiguous = params.keep_ambiguous_long ? '' : "${label}.ambiguous.ids"
    def keep_amb       = params.keep_ambiguous_long ? '--keep-ambiguous' : ''
    def cond_arg       = conditions.name != 'NO_FILE_COND' ? "--conditions ${conditions}" : ''
    // no genome: an empty host PAF, so every contaminant placement classifies as 'contaminant'
    def host_aln       = genome.name == 'NO_FILE_GENOME' ? ': > vs_host.paf'
                       : "minimap2 -x splice -t ${task.cpus} --secondary=no ${genome} ${seqs} > vs_host.paf 2> vs_host.log"
    """
    ${host_aln}
    minimap2 -x splice -t ${task.cpus} --secondary=no ${contaminants} ${seqs} > vs_contam.paf 2> vs_contam.log

    classify_host_contaminant.py vs_host.paf vs_contam.paf ${label} ${params.contaminant_margin}

    { cat ${label}.contaminant.ids ${drop_ambiguous} 2>/dev/null || true; } | sort -u > ${label}.removed.ids
    if [ -s ${label}.removed.ids ]; then
        seqkit grep -v -f ${label}.removed.ids ${seqs} > ${label}.hostonly.fasta
        seqkit grep    -f ${label}.removed.ids ${seqs} > ${label}.contaminant_removed.fasta
    else
        cp ${seqs} ${label}.hostonly.fasta
        : > ${label}.contaminant_removed.fasta
    fi

    N_IN=\$(grep -c '^>' ${seqs} || true)
    N_H=\$(wc -l < ${label}.host.ids); N_S=\$(wc -l < ${label}.contaminant.ids); N_A=\$(wc -l < ${label}.ambiguous.ids)
    {
      printf 'class\\tn_sequences\\n'
      printf 'host\\t%s\\n'            "\$N_H"
      printf 'no_alignment\\t%s\\n'    "\$(( N_IN - N_H - N_S - N_A ))"
      printf 'contaminant\\t%s\\n'     "\$N_S"
      printf 'ambiguous\\t%s\\n'       "\$N_A"
      printf 'kept\\t%s\\n'            "\$(grep -c '^>' ${label}.hostonly.fasta || true)"
      printf 'removed\\t%s\\n'         "\$(grep -c '^>' ${label}.contaminant_removed.fasta || true)"
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
    def label = meta.id
    """
    touch ${label}.hostonly.fasta ${label}.contaminant_removed.fasta ${label}.classification.tsv
    printf "class\\tn_sequences\\nhost\\t100\\n" > ${label}.longread_decontam.tsv
    printf "sample\\tcondition\\tisoforms_in\\trrna_removed\\thost\\tno_alignment\\tcontaminant\\tambiguous\\tcontaminant_removed\\tkept\\tpct_removed\\n" > longread_decontam_per_sample.tsv
    printf "LR1\\tcond\\t1000\\t10\\t800\\t100\\t80\\t10\\t90\\t900\\t9.09\\n" >> longread_decontam_per_sample.tsv
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}
