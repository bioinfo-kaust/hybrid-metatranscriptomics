/*
 * Differential expression with DESeq2 on tximport-aggregated salmon counts. The design formula,
 * the reference level and the contrast list are parameters, so the same code serves a
 * 3-treatment UV experiment and any other factorial design.
 */
process DESEQ2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:06e7afe89f992324'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:f4860cc0c268f684' }"

    input:
    path quants, stageAs: 'quants/*'
    tuple val(meta), path(tx2gene)
    path samplesheet
    path contrasts
    path expressed_tx

    output:
    tuple val(meta), path("dge_${meta.id}"), emit: results
    path "dge_${meta.id}/DE_summary.csv"   , emit: summary
    path "versions.yml"                    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label        = meta.id
    def contrast_arg = contrasts.name != 'NO_FILE_CONTRASTS' ? "--contrasts ${contrasts}" : ''
    // the DE-eligible set: genes with at least one transcript the experiment saw above the TPM
    // rule. The reference itself keeps the ORF+homology rescues so they stay in the index and in
    // the reported transcriptome; they just do not get tested. See filter_by_expression.py.
    def expr_tx_arg  = expressed_tx.name != 'NO_FILE_EXPRTX' ? "--expressed-tx ${expressed_tx}" : ''
    def ref_level    = params.reference_level ? "--reference-level ${params.reference_level}" : ''
    def rlibs        = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    deseq2.R --quant-dir quants --tx2gene ${tx2gene} --samplesheet ${samplesheet} \\
             --condition-column '${params.condition_column}' ${ref_level} \\
             --design '${params.design}' ${contrast_arg} \\
             --min-count ${params.de_min_count} --min-samples ${params.de_min_samples} \\
             --alpha ${params.de_alpha} --lfc ${params.de_lfc} --shrink ${params.de_shrink} \\
             ${expr_tx_arg} \\
             --outdir dge_${label} 2>&1 | tee deseq2.log
    mkdir -p dge_${label} && cp deseq2.log dge_${label}/
    # Ship the gene map with the results. Corset assigns cluster ids per run — only 14.6% of
    # transcripts keep the same Cluster-N.M label when salmon is re-run, though 99.91% of clusters
    # are identical as transcript sets — so a corset DE table is uninterpretable, and cannot be
    # joined to another run's, without the map that produced it. The locus map has the same
    # problem for a different reason: LOCUS_* ids are positional and shift when the reference does.
    cp ${tx2gene} dge_${label}/tx2gene.used.tsv
    emit_versions.sh "${task.process}" Rscript
    """

    stub:
    """
    mkdir -p dge_${meta.id} && touch dge_${meta.id}/DE_summary.csv
    emit_versions.sh "${task.process}" Rscript
    """
}
