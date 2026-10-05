/*
 * One per-transcript annotation table: ORFs, DIAMOND hits per database, eggNOG and InterPro.
 * InterProScan and eggNOG are optional; without them the table loses those columns only.
 */
process INTEGRATE_ANNOTATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(tx2gene), path(pep)
    path diamond_hits, stageAs: 'diamond/*'
    path emapper
    path interpro

    output:
    tuple val(meta), path("annotation.${meta.id}.tsv"), emit: table
    path "annotation.${meta.id}.summary.txt"          , emit: summary
    path "versions.yml"                               , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label = meta.id
    def egg   = emapper.name  != 'NO_FILE'     ? "--emapper ${emapper}"   : ''
    def ipr   = interpro.name != 'NO_FILE_IPR' ? "--interpro ${interpro}" : ''
    def prim  = params.primary_protein_db ? "--primary-db ${params.primary_protein_db}" : ''
    """
    integrate_annotation.py \\
        --tx2gene ${tx2gene} --pep ${pep} \\
        --diamond-dir diamond ${egg} ${ipr} ${prim} \\
        --interest-regex '${params.interest_regex}' \\
        --interest-label '${params.interest_label}' \\
        --out annotation.${label}.tsv \\
        --summary annotation.${label}.summary.txt
    cat annotation.${label}.summary.txt
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch annotation.${meta.id}.tsv annotation.${meta.id}.summary.txt
    emit_versions.sh "${task.process}" python3
    """
}
