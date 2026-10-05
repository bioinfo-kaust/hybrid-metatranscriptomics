/* Per-sample salmon mapping rates; fails the run when a library maps below --min_mapping_rate. */
process QUANT_SUMMARY {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    path quants, stageAs: 'quants/*'

    output:
    path 'mapping_rates.tsv', emit: table
    path "versions.yml"     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    quant_summary.py --dir quants --out mapping_rates.tsv --min-mapping-rate ${params.min_mapping_rate}
    cat mapping_rates.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch mapping_rates.tsv
    emit_versions.sh "${task.process}" python3
    """
}
