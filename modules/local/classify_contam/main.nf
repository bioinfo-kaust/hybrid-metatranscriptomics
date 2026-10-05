/*
 * Residual audit of an assembled set. After read-level removal this is a QC gate, not the
 * separation step: it reports how much contaminant signal survived, and drops what did.
 * keep_class controls the keep rule ('host,ambiguous,unmapped' in this pipeline).
 */
process CLASSIFY_CONTAM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:49f77f39b14c7e4f'
        : 'community.wave.seqera.io/library/minimap2_samtools_seqkit_python_gawk:9b379af53c51bfba' }"

    input:
    tuple val(meta), path(fasta)
    path genome
    path contaminants
    val keep_class

    output:
    tuple val(meta), path("${meta.id}.host.fasta"), emit: fasta
    path "${meta.id}.hostsym_counts.txt"          , emit: counts
    path "${meta.id}.classification.tsv"          , emit: classification
    path "${meta.id}.host.ids"                    , emit: host_ids
    path "versions.yml"                           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label    = meta.id
    def host_aln = genome.name == 'NO_FILE_GENOME' ? ': > vs_host.paf'
                 : "minimap2 -x splice -t ${task.cpus} --secondary=no ${genome} ${fasta} > vs_host.paf 2> vs_host.log"
    """
    ${host_aln}
    minimap2 -x splice -t ${task.cpus} --secondary=no ${contaminants} ${fasta} > vs_contam.paf 2> vs_contam.log
    classify_host_contaminant.py vs_host.paf vs_contam.paf ${label} ${params.contaminant_margin}

    # ids up to the first whitespace; seqkit reads the input gzipped or not
    seqkit seq -n -i ${fasta} | sort -u > all.ids

    # A sequence that aligns to NEITHER reference never appears in either PAF, so the classifier
    # never lists it in any class — it is not "unmapped: 0", it is absent. Derive that class from
    # the FASTA instead, or the keep rule deletes it: 19,166 of 129,450 Trinity contigs (14.8%)
    # and 2,804 of 79,140 reference-free long-read sequences on this dataset, which is exactly the
    # genome-absent material the de novo route exists to recover.
    { cat ${label}.host.ids ${label}.contaminant.ids ${label}.ambiguous.ids 2>/dev/null || true; } \\
      | sort -u > classified.ids
    comm -23 all.ids classified.ids > ${label}.unmapped.ids

    : > keep.ids
    for c in \$(echo "${keep_class}" | tr ',' ' '); do
        # `[ -s f ] && cat f` as the last command in the body exits non-zero on an empty class
        # file, which under `bash -e` aborts the task
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
      echo "  kept (${keep_class}): \$(grep -c '^>' ${label}.host.fasta || true)"
      echo "  dropped: \$(( \$(wc -l < all.ids) - \$(grep -c '^>' ${label}.host.fasta || true) ))"
    } | tee ${label}.hostsym_counts.txt
    emit_versions.sh "${task.process}" minimap2 seqkit
    """

    stub:
    def label = meta.id
    """
    touch ${label}.host.fasta ${label}.hostsym_counts.txt ${label}.classification.tsv ${label}.host.ids
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}
