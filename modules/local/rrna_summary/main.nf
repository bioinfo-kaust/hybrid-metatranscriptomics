/*
 * Per-sample rRNA percentage, parsed from bbduk's own accounting (BBMAP_BBDUK's log).
 * --paired: bbduk counts reads, every other per-sample table here counts pairs.
 */
process RRNA_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(bbduk_log)

    output:
    path "${meta.id}.rrna.tsv", emit: summary
    path "versions.yml"       , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    rrna_summary.py ${meta.id} ${bbduk_log} ${meta.condition} --paired > ${meta.id}.rrna.tsv
    cat ${meta.id}.rrna.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    printf "sample\\tcondition\\treads_in\\trrna_removed\\treads_kept\\trrna_pct\\tunit\\n${meta.id}\\t${meta.condition}\\t1000\\t100\\t900\\t10.00\\tpairs\\n" > ${meta.id}.rrna.tsv
    emit_versions.sh "${task.process}" python3
    """
}
