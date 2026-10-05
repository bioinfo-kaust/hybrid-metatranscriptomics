/* TransDecoder marks stop codons with '*', which several downstream tools reject. */
process CLEAN_PEP {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(pep)

    output:
    tuple val(meta), path("${meta.id}.clean.pep"), emit: pep
    path "versions.yml"                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sed '/^>/!s/\\*//g' ${pep} > ${meta.id}.clean.pep
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch ${meta.id}.clean.pep
    emit_versions.sh "${task.process}"
    """
}
