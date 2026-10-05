/*
 * The interactive dashboard: one self-contained HTML file built from the run's own outputs.
 * It carries results and methods together — every stage card shows the arguments it actually ran
 * with (read from the run manifest) beside the numbers it produced.
 */
process DASHBOARD {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    path inputs, stageAs: 'in/*'
    path manifest
    val title

    output:
    path 'dashboard/index.html', emit: html
    path "versions.yml"        , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def man       = manifest.name != 'NO_FILE_MANIFEST' ? "--manifest ${manifest}" : ''
    def safe_title = title.toString().replace("'", "")
    """
    build_dashboard.py --dir in ${man} --out dashboard/index.html --title '${safe_title}'
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    mkdir -p dashboard && touch dashboard/index.html
    emit_versions.sh "${task.process}" python3
    """
}
