/*
 * Reporting — aggregated QC tables, MultiQC, and the figure set.
 */

process AGGREGATE_TSV {
    label 'process_single'
    tag "$name"

    input:
    tuple val(name), path(tables, stageAs: 'tab/*')

    output:
    tuple val(name), path("${name}.tsv"), emit: table
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    # keep the header of the first file only
    first=\$(ls tab/* | head -1)
    head -1 "\$first" > ${name}.tsv
    for f in tab/*; do tail -n +2 "\$f" >> ${name}.tsv; done
    column -t ${name}.tsv || cat ${name}.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch ${name}.tsv
    emit_versions.sh "${task.process}" python3
    """
}

/* Per-sample decontamination table -> a MultiQC custom-content bargraph. */
process DECONTAM_MQC {
    label 'process_single'

    input:
    path table

    output:
    path 'decontam_mqc.tsv', emit: mqc
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    {
      echo "# id: 'contaminant_removal'"
      echo "# section_name: 'Contaminant read removal'"
      echo "# description: 'Read pairs whose best genome-wide placement was a contaminant scaffold (removed) versus pairs retained as host (host placement + reads matching neither reference).'"
      echo "# plot_type: 'bargraph'"
      echo "# pconfig:"
      echo "#     id: 'contaminant_bargraph'"
      echo "#     title: 'Contaminant vs host read pairs per sample'"
      echo "#     ylab: 'read pairs'"
      printf 'Sample\\tContaminant_removed\\tHost_kept\\n'
      awk -F'\\t' 'NR>1{printf "%s\\t%s\\t%s\\n", \$1, \$4, \$5}' ${table}
    } > decontam_mqc.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch decontam_mqc.tsv
    emit_versions.sh "${task.process}" python3
    """
}

process MULTIQC {
    label 'process_medium'

    input:
    path '*'
    path config

    output:
    path '*multiqc_report.html'          , emit: report
    path '*_data'                        , emit: data
    path '*_plots' , optional: true      , emit: plots
    path "versions.yml", emit: versions, topic: versions

    script:
    def cfg = config.name != 'NO_FILE_MQC' ? "-c ${config}" : ''
    """
    multiqc -f -o . -n multiqc_report.html ${cfg} .
    emit_versions.sh "${task.process}" multiqc
    """

    stub:
    """
    # MultiQC names the data directory after the report: multiqc_report.html -> multiqc_report_data
    mkdir -p multiqc_report_data && touch multiqc_report.html
    emit_versions.sh "${task.process}" multiqc
    """
}

/*
 * Figure set. Every panel is computed from this run's own outputs — nothing is hard-coded,
 * and a panel whose input is missing is skipped rather than failing the run.
 */
process FIGURES {
    label 'process_medium'
    label 'error_ignore'
    tag "$label"

    input:
    tuple val(label), path(dge_dir), path(enr_dir)
    tuple val(label3), path(annotation)
    path qc_tables, stageAs: 'qc/*'
    path busco    , stageAs: 'busco/*'
    path ladder

    output:
    path "figures_${label}", emit: figures
    path "versions.yml", emit: versions, topic: versions

    script:
    def rlibs = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    mkdir -p figures_${label}
    figures.R --dge-dir ${dge_dir} --enrichment-dir ${enr_dir} --annotation ${annotation} \\
              --qc-dir qc --busco-dir busco --ladder ${ladder} \\
              --condition-column '${params.condition_column}' \\
              --padj ${params.de_alpha} --lfc ${params.de_lfc} \\
              --outdir figures_${label} 2>&1 | tee figures_${label}/figures.log
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    mkdir -p figures_${label} && touch figures_${label}/F1.png
    emit_versions.sh "${task.process}" python3
    """
}

/*
 * The interactive dashboard: one self-contained HTML file built from the run's own outputs.
 * It carries results and methods together — every stage card shows the arguments it actually ran
 * with (read from the run manifest) beside the numbers it produced.
 */
process DASHBOARD {
    label 'process_low'

    input:
    path inputs, stageAs: 'in/*'
    path manifest

    output:
    path 'dashboard/index.html', emit: html
    path "versions.yml", emit: versions, topic: versions

    script:
    def man   = manifest.name != 'NO_FILE_MANIFEST' ? "--manifest ${manifest}" : ''
    def title = (params.dashboard_title ?: "Transcriptome run — ${params.dge_reference}").replace("'", "")
    """
    build_dashboard.py --dir in ${man} --out dashboard/index.html --title '${title}'
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    mkdir -p dashboard && touch dashboard/index.html
    emit_versions.sh "${task.process}" python3
    """
}

/* A single human-readable summary of what the run produced, with the numbers to check. */
process RUN_SUMMARY {
    label 'process_single'

    input:
    path tables, stageAs: 'in/*'

    output:
    path 'run_summary.txt', emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    run_summary.py --dir in --out run_summary.txt --primary-map ${Inputs.primaryMap(params)}
    cat run_summary.txt
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch run_summary.txt
    emit_versions.sh "${task.process}" python3
    """
}
