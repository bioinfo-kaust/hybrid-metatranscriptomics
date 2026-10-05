/*
 * Differential isoform usage. Structure comes from the long-read collapsed GTF, abundances from
 * the short-read quantification — so the GTF must be restricted to the long-read models that
 * survive into THIS reference set, or importRdata's Jaccard check fails.
 */
process ISOFORM_SWITCHING {
    tag "${meta.id}"
    label 'process_high'
    label 'error_ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bioconductor-isoformswitchanalyzer_gawk:c5fe3261a639c47b'
        : 'community.wave.seqera.io/library/bioconductor-isoformswitchanalyzer_gawk:27d6b2549f174448' }"

    input:
    path quants, stageAs: 'quants/*'
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)
    path samplesheet

    output:
    tuple val(meta), path("isoforms_${meta.id}"), emit: results, optional: true
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def label = meta.id
    def rlibs = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    { grep '^>' ${fasta} || true; } | sed 's/^>//;s/[ |].*//' | sort -u > set_ids.txt
    awk -v ids=set_ids.txt '
      BEGIN{ while ((getline l < ids) > 0) keep[l]=1 }
      match(\$0, /transcript_id "[^"]+"/) {
        t = substr(\$0, RSTART+15, RLENGTH-16); if (t in keep) print
      }' ${gtf} > subset.gtf
    echo "[iso] GTF restricted to ${label}: \$(awk -F'\\t' '\$3=="transcript"' subset.gtf | wc -l) of \$(awk -F'\\t' '\$3=="transcript"' ${gtf} | wc -l) models"

    mkdir -p isoforms_${label}
    isoform_switching.R --quant-dir quants --fasta ${fasta} --gtf subset.gtf \\
        --samplesheet ${samplesheet} --condition-column '${params.condition_column}' \\
        --outdir isoforms_${label} 2>&1 | tee isoforms_${label}/isoform_switching.log
    emit_versions.sh "${task.process}" Rscript
    """

    stub:
    """
    mkdir -p isoforms_${meta.id} && touch isoforms_${meta.id}/switch_summary.csv
    emit_versions.sh "${task.process}" Rscript
    """
}
