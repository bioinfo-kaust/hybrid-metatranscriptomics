/*
 * Figure set. Every panel is computed from this run's own outputs — nothing is hard-coded,
 * and a panel whose input is missing is skipped rather than failing the run.
 */
process FIGURES {
    tag "${meta.id}"
    label 'process_medium'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:06e7afe89f992324'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:f4860cc0c268f684' }"

    input:
    tuple val(meta), path(dge_dir), path(enr_dir)
    tuple val(meta2), path(annotation)
    path qc_tables, stageAs: 'qc/*'
    path busco    , stageAs: 'busco/*'
    path ladder

    output:
    path "figures_${meta.id}", emit: figures
    path "versions.yml"      , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label = meta.id
    def rlibs = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    mkdir -p figures_${label}
    figures.R --dge-dir ${dge_dir} --enrichment-dir ${enr_dir} --annotation ${annotation} \\
              --qc-dir qc --busco-dir busco --ladder ${ladder} \\
              --condition-column '${params.condition_column}' \\
              --padj ${params.de_alpha} --lfc ${params.de_lfc} \\
              --outdir figures_${label} 2>&1 | tee figures_${label}/figures.log
    emit_versions.sh "${task.process}" Rscript
    """

    stub:
    """
    mkdir -p figures_${meta.id} && touch figures_${meta.id}/F1.png
    emit_versions.sh "${task.process}" Rscript
    """
}
