/* Concatenate per-sample tables into one, keeping the first file's header. */
process AGGREGATE_TSV {
    tag "${name}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(name), path(tables, stageAs: 'tab/*')

    output:
    tuple val(name), path("${name}.tsv"), emit: table
    path "versions.yml"                 , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    files=( tab/* )
    head -1 "\${files[0]}" > ${name}.tsv
    for f in "\${files[@]}"; do tail -n +2 "\$f" >> ${name}.tsv; done
    column -t ${name}.tsv 2>/dev/null || cat ${name}.tsv
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch ${name}.tsv
    emit_versions.sh "${task.process}"
    """
}
