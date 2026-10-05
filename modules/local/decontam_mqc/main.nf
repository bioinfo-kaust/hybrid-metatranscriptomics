/* Per-sample decontamination table -> a MultiQC custom-content bargraph. */
process DECONTAM_MQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    path table

    output:
    path 'decontam_mqc.tsv', emit: mqc
    path "versions.yml"    , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

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
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    touch decontam_mqc.tsv
    emit_versions.sh "${task.process}"
    """
}
