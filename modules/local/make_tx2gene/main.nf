/*
 * The string-stripped transcript -> gene map: per-source rules ('prefix:isoform_suffix_regex',
 * first match wins) applied to the transcript ids.
 */
process MAKE_TX2GENE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("tx2gene.${meta.id}.tsv"), emit: tx2gene
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    make_tx2gene.py --fasta ${fasta} --rules '${params.tx2gene_rules}' --out tx2gene.${meta.id}.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    printf "t1\\tg1\\n" > tx2gene.${meta.id}.tsv
    emit_versions.sh "${task.process}" python3
    """
}
