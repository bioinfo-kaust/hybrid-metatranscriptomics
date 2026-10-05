/* A single human-readable summary of what the run produced, with the numbers to check. */
process RUN_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    path tables, stageAs: 'in/*'
    val primary_map

    output:
    path 'run_summary.txt', emit: summary
    path "versions.yml"   , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    run_summary.py --dir in --out run_summary.txt --primary-map ${primary_map}
    cat run_summary.txt
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch run_summary.txt
    emit_versions.sh "${task.process}" python3
    """
}
