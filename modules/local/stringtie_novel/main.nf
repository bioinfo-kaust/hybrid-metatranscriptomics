/*
 * Keep only the merged short-read loci the long-read models missed, and give them sequence.
 * gffcompare class codes u/x/i/y = intergenic / antisense / intronic / contained-in-intron;
 * everything else is a variant of a locus the long reads already describe. Without long-read
 * models to compare against (annotated = NO_FILE_GFFCMP) every merged transcript is taken.
 */
process STRINGTIE_NOVEL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/gffread_seqkit_gawk:97e9b64ec67d60ea'
        : 'community.wave.seqera.io/library/gffread_seqkit_gawk:7d0bb4806e001c3b' }"

    input:
    tuple val(meta), path(merged_gtf)
    path annotated_gtf
    path genome

    output:
    tuple val(meta), path('stringtie_novel.fasta'), emit: fasta
    path 'stringtie_novel_ids.txt'                , emit: novel_ids
    path "versions.yml"                           , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def codes   = params.stringtie_class_codes
    def has_ref = annotated_gtf.name != 'NO_FILE_GFFCMP'
    """
    if ${has_ref}; then
        # POSIX awk: match() takes two arguments, and the 3-argument gawk form used to fail
        # silently, leaving an empty StringTie contribution that looked like a real result.
        # RSTART/RLENGTH carry the same information portably.
        awk -v codes="${codes}" '
            \$3 == "transcript" {
                cc = ""; tid = ""
                if (match(\$0, /class_code "[^"]+"/)) cc  = substr(\$0, RSTART+12, RLENGTH-13)
                if (match(\$0, /transcript_id "[^"]+"/)) tid = substr(\$0, RSTART+15, RLENGTH-16)
                if (cc != "" && tid != "" && cc ~ "^(" codes ")\$") print tid
            }' ${annotated_gtf} | sort -u > stringtie_novel_ids.txt
        # select records by PARSED transcript_id: -Fw matched the id as a substring anywhere on
        # the line, so a gene_id sharing the prefix pulled in records that were never selected.
        awk 'NR==FNR{keep[\$1]=1; next}
             { if (match(\$0, /transcript_id "[^"]+"/)) {
                 t = substr(\$0, RSTART+15, RLENGTH-16); if (t in keep) print }
             }' stringtie_novel_ids.txt ${annotated_gtf} > novel.gtf
    else
        # no long-read models to compare against: every merged transcript is the base set
        { grep -oE 'transcript_id "[^"]+"' ${merged_gtf} || true; } | sed 's/transcript_id "//;s/"//' \\
          | sort -u > stringtie_novel_ids.txt
        cp ${merged_gtf} novel.gtf
    fi
    echo "short-read transcripts taken forward: \$(wc -l < stringtie_novel_ids.txt)"
    gffread -w novel.raw.fasta -g ${genome} novel.gtf
    awk '/^>/{sub(/^>/,">${params.prefix_stringtie}")}1' novel.raw.fasta \\
      | seqkit seq -m ${params.min_transcript_len} - > stringtie_novel.fasta
    echo "novel >=${params.min_transcript_len} bp: \$(grep -c '^>' stringtie_novel.fasta || true)"
    emit_versions.sh "${task.process}" gffread seqkit
    """

    stub:
    """
    touch stringtie_novel.fasta stringtie_novel_ids.txt
    emit_versions.sh "${task.process}" gffread seqkit
    """
}
