/*
 * Restrict a GTF to the transcripts named in an id list. Used to score StringTie loci against the
 * long-read models that SURVIVED the contaminant audit: comparing against the raw collapse output
 * meant a StringTie locus matching a model later deleted as contaminant was judged "not novel"
 * and disappeared from every rung.
 */
process FILTER_GTF_BY_IDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(gtf)
    path ids

    output:
    tuple val(meta), path("${meta.id}.audited.gtf"), emit: gtf
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # select records by PARSED transcript_id, not by substring
    awk 'NR==FNR{keep[\$1]=1; next}
         { if (match(\$0, /transcript_id "[^"]+"/)) {
             t = substr(\$0, RSTART+15, RLENGTH-16); if (t in keep) print }
         }' ${ids} ${gtf} > ${meta.id}.audited.gtf
    echo "reference models: \$(grep -c 'transcript_id' ${gtf} || true) -> \$(grep -c 'transcript_id' ${meta.id}.audited.gtf || true) after the audit"
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch ${meta.id}.audited.gtf
    emit_versions.sh "${task.process}"
    """
}
