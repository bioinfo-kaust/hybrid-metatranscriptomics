/* Restrict a tx2gene map to the transcripts that survive into a deduplicated set. */
process SUBSET_TX2GENE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(fasta)
    path tx2gene

    output:
    tuple val(meta), path("tx2gene.${meta.id}.tsv"), emit: tx2gene
    path "versions.yml"                            , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    { grep '^>' ${fasta} || true; } | sed 's/^>//; s/[ |].*//' | sort -u > ids.txt
    awk 'NR==FNR{k[\$1];next} (\$1 in k)' ids.txt ${tx2gene} > tx2gene.${meta.id}.tsv
    echo "[${meta.id}] transcripts \$(wc -l < ids.txt), genes \$(cut -f2 tx2gene.${meta.id}.tsv | sort -u | wc -l)"
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    printf "t1\\tg1\\n" > tx2gene.${meta.id}.tsv
    emit_versions.sh "${task.process}" python3
    """
}
