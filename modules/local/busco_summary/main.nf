/*
 * The one-line completeness string from a BUSCO short summary, as <label>.busco.txt — the form the
 * ladder table, the figures and the dashboard read.
 *
 * It also carries BUSCO's full_table.tsv out, renamed to <label>.full_table.tsv. That table says
 * which ortholog each sequence carries, and it is the label-free truth an external gene-map
 * identity evaluation scores the tx2gene maps against. It cannot be published straight from
 * BUSCO_BUSCO: that process declares the whole `*-busco` directory as an output as well as the
 * nested table, and Nextflow then treats the directory as the unit of publication — a publishDir
 * pattern matching the nested file publishes nothing (verified), and publishing the directory
 * instead would cost 123 MB per rung for a 512 KB table. Copying it here, where the ids are
 * already being translated for downstream consumers, is both cheaper and the same kind of work.
 */
process BUSCO_SUMMARY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_gawk:b4b1f7b00b38c493'
        : 'community.wave.seqera.io/library/python_gawk:d0e9449be7e08c97' }"

    input:
    tuple val(meta), path(summary, stageAs: 'summary/*'), path(full_table, stageAs: 'ft/*')

    output:
    tuple val(meta), path("${meta.id}.busco.txt")      , emit: summary
    // absent when BUSCO failed for this rung (errorStrategy ignore), hence optional
    tuple val(meta), path("${meta.id}.full_table.tsv") , emit: full_table, optional: true
    path "versions.yml"                                , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    { grep -hoE "C:[0-9.]+%\\[S:[0-9.]+%,D:[0-9.]+%\\],F:[0-9.]+%,M:[0-9.]+%" summary/*.txt 2>/dev/null || true; } \\
      | sed -n 1p > ${meta.id}.busco.txt
    [ -s ${meta.id}.busco.txt ] || echo NA >> ${meta.id}.busco.txt
    echo "[${meta.id}] \$(cat ${meta.id}.busco.txt)"

    # ft/ holds full_table.tsv when BUSCO produced one; no table is not an error, it just means
    # this rung has no ortholog assignment to report
    for f in ft/full_table.tsv ft/*.tsv; do
        if [ -s "\$f" ]; then cp "\$f" ${meta.id}.full_table.tsv; break; fi
    done
    if [ -s ${meta.id}.full_table.tsv ]; then
        echo "[${meta.id}] full table: \$(grep -vc '^#' ${meta.id}.full_table.tsv) rows"
    else
        echo "[${meta.id}] no BUSCO full table staged"
    fi
    emit_versions.sh "${task.process}"
    """

    stub:
    """
    echo "C:70.0%[S:50.0%,D:20.0%],F:10.0%,M:20.0%" > ${meta.id}.busco.txt
    printf '# Busco id\\tStatus\\tSequence\\n1at6073\\tComplete\\tt1\\n' > ${meta.id}.full_table.tsv
    emit_versions.sh "${task.process}"
    """
}
