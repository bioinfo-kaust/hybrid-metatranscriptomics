/*
 * Count matrices for every gene map, and for the transcripts underneath them.
 *
 * The run already publishes per-sample salmon quant.sf and, from DESeq2, normalised and VST
 * counts -- but those are transformed, and they cover only the genes that passed the presence
 * pre-filter. Nothing published the raw gene-level matrix, so the gene FASTA shipped without the
 * counts belonging to it. This writes that matrix, plus the length offset DESeq2 pairs it with,
 * a length-corrected variant for tools that cannot take an offset, TPM, and the same three
 * quantities at transcript level.
 *
 * It re-runs the aggregation rather than re-quantifying, and it is handed each map's dds.rds so it
 * can prove the result IS the matrix DESeq2 modelled -- see bin/count_matrix.R. One task for all
 * maps, because the transcript-level matrices are map-independent and should exist once.
 */
process COUNT_MATRIX {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:06e7afe89f992324'
        : 'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_bioconductor-apeglm_r-ashr_pruned:f4860cc0c268f684' }"

    input:
    path quants, stageAs: 'quants/*'
    path tx2gene_files
    path gene_fastas
    path transcript_fasta
    path dge_dirs
    path samplesheet

    output:
    path "gene_level"   , emit: matrices
    path "versions.yml" , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def t2g   = tx2gene_files.collect { it.name }.join(',')
    def fas   = gene_fastas ? gene_fastas.collect { it.name }.join(',') : ''
    // the dge_* directories, each holding the dds.rds that verifies its own map's matrix
    def dges  = dge_dirs ? dge_dirs.collect { it.name }.join(',') : ''
    def rlibs = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    count_matrix.R --quant-dir quants \\
                   --samplesheet ${samplesheet} \\
                   --tx2gene '${t2g}' \\
                   ${fas ? "--gene-fasta '${fas}'" : ''} \\
                   ${dges ? "--dge-dir '${dges}'" : ''} \\
                   --transcript-fasta ${transcript_fasta} \\
                   --primary-map '${params.dge_primary_map}' \\
                   --outdir gene_level 2>&1 | tee count_matrix.log
    cp count_matrix.log gene_level/
    emit_versions.sh "${task.process}" Rscript
    """

    stub:
    """
    mkdir -p gene_level && touch gene_level/checks.txt
    emit_versions.sh "${task.process}" Rscript
    """
}
