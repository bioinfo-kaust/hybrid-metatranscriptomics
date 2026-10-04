/*
 * Stage 4 — functional annotation of the reference transcriptome.
 *   TransDecoder (ORFs/peptides) -> DIAMOND vs one or more protein DBs + eggNOG-mapper
 *                                -> InterProScan (slow, optional)
 *                                -> one integrated per-transcript table
 */

process TRANSDECODER {
    label 'process_medium'
    tag "$label"

    input:
    tuple val(label), path(fasta)

    output:
    tuple val(label), path("${label}.pep"), emit: pep
    tuple val(label), path("${label}.cds"), emit: cds
    path "${label}.transdecoder.gff3"     , emit: gff3, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '--single_best_only'
    """
    # TransDecoder resumes from checkpoints purely on file presence, so a stale working dir
    # silently emits peptides keyed on a previous reference's ids. Nextflow gives each task a
    # clean directory, which removes that failure mode by construction.
    TransDecoder.LongOrfs -t ${fasta}
    TransDecoder.Predict -t ${fasta} ${args} --no_refine_starts 2> predict.log \\
      || TransDecoder.Predict -t ${fasta} ${args} 2>> predict.log

    cp *.transdecoder.pep  ${label}.pep
    cp *.transdecoder.cds  ${label}.cds
    cp *.transdecoder.gff3 ${label}.transdecoder.gff3 2>/dev/null || true
    echo "[${label}] ORFs: \$(grep -c '^>' ${label}.pep)"
    emit_versions.sh "${task.process}" TransDecoder.LongOrfs
    """

    stub:
    """
    touch ${label}.pep ${label}.cds ${label}.transdecoder.gff3
    emit_versions.sh "${task.process}" TransDecoder.LongOrfs
    """
}

/* TransDecoder marks stop codons with '*', which several downstream tools reject. */
process CLEAN_PEP {
    label 'process_single'
    tag "$label"

    input:
    tuple val(label), path(pep)

    output:
    tuple val(label), path("${label}.clean.pep"), emit: pep
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    sed '/^>/!s/\\*//g' ${pep} > ${label}.clean.pep
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch ${label}.clean.pep
    emit_versions.sh "${task.process}" python3
    """
}

process DIAMOND_BLASTP {
    label 'process_high'
    tag "${label} vs ${db_name}"

    input:
    tuple val(label), path(pep), val(db_name), path(db)

    output:
    tuple val(db_name), path("diamond_${db_name}.tsv"), emit: hits
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '-k 1 -e 1e-5 --more-sensitive'
    """
    diamond blastp -q ${pep} -d ${db} -o diamond_${db_name}.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \\
        ${args} -p ${task.cpus} --quiet
    echo "[${db_name}] hits: \$(wc -l < diamond_${db_name}.tsv)"
    emit_versions.sh "${task.process}" diamond
    """

    stub:
    """
    touch diamond_${db_name}.tsv
    emit_versions.sh "${task.process}" diamond
    """
}

process EGGNOG_MAPPER {
    label 'process_high'
    tag "$label"

    input:
    tuple val(label), path(pep)
    val data_dir

    output:
    path 'eggnog.emapper.annotations', emit: annotations
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '-m diamond --override'
    """
    emapper.py -i ${pep} --itype proteins -o eggnog --output_dir . \\
        --data_dir ${data_dir} --cpu ${task.cpus} ${args} > emapper.log 2>&1
    echo "[eggnog] rows: \$(grep -vc '^#' eggnog.emapper.annotations)"
    emit_versions.sh "${task.process}" emapper.py
    """

    stub:
    """
    touch eggnog.emapper.annotations
    emit_versions.sh "${task.process}" emapper.py
    """
}

process INTERPROSCAN {
    label 'process_high'
    label 'process_long'
    tag "$label"

    input:
    tuple val(label), path(pep)

    output:
    path 'interpro.tsv', emit: tsv
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '-dp -goterms'
    """
    interproscan.sh -i ${pep} -f tsv -o interpro.tsv \\
        -appl ${params.interproscan_appl} -cpu ${task.cpus} ${args} \\
        > /dev/null 2> interproscan.log \\
      || { echo "InterProScan failed"; tail -20 interproscan.log; exit 1; }
    echo "[interpro] rows: \$(wc -l < interpro.tsv)"
    emit_versions.sh "${task.process}" interproscan.sh
    """

    stub:
    """
    touch interpro.tsv
    emit_versions.sh "${task.process}" interproscan.sh
    """
}

process INTEGRATE_ANNOTATION {
    label 'process_medium'
    tag "$label"

    input:
    tuple val(label), path(tx2gene), path(pep)
    path diamond_hits, stageAs: 'diamond/*'
    path emapper
    path interpro

    output:
    tuple val(label), path("annotation.${label}.tsv"), emit: table
    path "annotation.${label}.summary.txt"           , emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    def egg = emapper.name != 'NO_FILE'      ? "--emapper ${emapper}"   : ''
    def ipr = interpro.name != 'NO_FILE_IPR' ? "--interpro ${interpro}" : ''
    """
    integrate_annotation.py \\
        --tx2gene ${tx2gene} --pep ${pep} \\
        --diamond-dir diamond ${egg} ${ipr} \\
        ${params.primary_protein_db ? "--primary-db ${params.primary_protein_db}" : ''} \\
        --interest-regex '${params.interest_regex}' \\
        --interest-label '${params.interest_label}' \\
        --out annotation.${label}.tsv \\
        --summary annotation.${label}.summary.txt
    cat annotation.${label}.summary.txt
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch annotation.${label}.tsv annotation.${label}.summary.txt
    emit_versions.sh "${task.process}" python3
    """
}
