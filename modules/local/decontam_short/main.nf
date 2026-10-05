/*
 * Short reads: one best placement per read on the combined host + contaminant index; a pair is
 * removed when either mate's primary alignment is on a contaminant scaffold.
 */
process DECONTAM_SHORT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:49f77f39b14c7e4f'
        : 'community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:9b379af53c51bfba' }"

    input:
    tuple val(meta), path(reads)
    path index
    val contaminant_prefix

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.host.fq.gz"), emit: reads
    path "${meta.id}.decontam.tsv"                       , emit: summary
    path "${meta.id}.minimap2.log"                       , emit: log
    path "versions.yml"                                  , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # -F 0x900 keeps primary, non-supplementary records: one best placement per read
    minimap2 -ax sr -t ${task.cpus} ${index} ${reads[0]} ${reads[1]} 2> ${meta.id}.minimap2.log \\
      | samtools view -F 0x900 - \\
      | awk -F'\\t' '\$3 ~ /^${contaminant_prefix}/ {print \$1}' \\
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
    echo "" | gzip > ${meta.id}_R1.host.fq.gz
    echo "" | gzip > ${meta.id}_R2.host.fq.gz
    touch ${meta.id}.minimap2.log
    printf "sample\\tcondition\\tpairs_in\\tcontaminant_removed\\tpairs_kept\\tpct_removed\\n${meta.id}\\t${meta.condition}\\t1000\\t130\\t870\\t13.00\\n" > ${meta.id}.decontam.tsv
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}
