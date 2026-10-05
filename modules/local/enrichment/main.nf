/* GO + KEGG over-representation and GO GSEA per contrast, against the pipeline's own annotation. */
process ENRICHMENT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:06e7afe89f992324'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:f4860cc0c268f684' }"

    input:
    tuple val(meta), path(dge_dir)
    tuple val(meta2), path(annotation)

    output:
    tuple val(meta), path("enrichment_${meta.id}"), emit: results
    path "versions.yml"                           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label = meta.id
    def rlibs = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    enrichment.R --annotation ${annotation} --dge-dir ${dge_dir} --outdir enrichment_${label} \\
                 --padj ${params.de_alpha} --lfc ${params.de_lfc} \\
                 --seed ${params.random_seed} 2>&1 | tee enrichment.log
    mkdir -p enrichment_${label} && cp enrichment.log enrichment_${label}/
    emit_versions.sh "${task.process}" Rscript
    """

    stub:
    """
    mkdir -p enrichment_${meta.id} && touch enrichment_${meta.id}/enrichment_summary.csv
    emit_versions.sh "${task.process}" Rscript
    """
}
