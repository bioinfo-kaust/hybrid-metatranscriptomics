/* The auditable ladder: one row per reference variant with transcripts, genes and BUSCO. */
process ROUTE_CONTRIBUTION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    path stats  , stageAs: 'stats/*'
    path tx2gene, stageAs: 't2g/*'
    path busco  , stageAs: 'busco/*'

    output:
    path 'route_contribution.tsv', emit: table
    path "versions.yml"          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    route_contribution.py --stats stats --tx2gene t2g --busco busco --out route_contribution.tsv
    column -t -s\$'\\t' route_contribution.tsv 2>/dev/null || cat route_contribution.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch route_contribution.tsv
    emit_versions.sh "${task.process}" python3
    """
}
