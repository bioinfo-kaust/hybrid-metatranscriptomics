/*
 * Stage 6 — differential expression, enrichment, differential isoform usage, and the
 * coordinate-based rebuild of gene identity.
 *
 * The design formula, the reference level and the contrast list are all parameters, so the
 * same code serves a 3-treatment UV experiment and any other factorial design.
 */

process DESEQ2 {
    label 'process_medium'
    tag "$label"

    input:
    path quants, stageAs: 'quants/*'
    tuple val(label), path(tx2gene)
    path samplesheet
    path contrasts
    path expressed_tx

    output:
    tuple val(label), path("dge_${label}"), emit: results
    path "dge_${label}/DE_summary.csv"    , emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    def contrast_arg = contrasts.name != 'NO_FILE_CONTRASTS' ? "--contrasts ${contrasts}" : ''
    // the DE-eligible set: genes with at least one transcript the experiment saw above the TPM
    // rule. The reference itself keeps the ORF+homology rescues so they stay in the index and in
    // the reported transcriptome; they just do not get tested. See filter_by_expression.py.
    def expr_tx_arg = expressed_tx.name != 'NO_FILE_EXPRTX' ? "--expressed-tx ${expressed_tx}" : ''
    def rlibs = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    deseq2.R --quant-dir quants --tx2gene ${tx2gene} --samplesheet ${samplesheet} \\
             --condition-column '${params.condition_column}' \\
             ${params.reference_level ? "--reference-level ${params.reference_level}" : ''} \\
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
    mkdir -p dge_${label} && touch dge_${label}/DE_summary.csv
    emit_versions.sh "${task.process}" Rscript
    """
}

process ENRICHMENT {
    label 'process_medium'
    tag "$label"

    input:
    tuple val(label), path(dge_dir)
    tuple val(label2), path(annotation)

    output:
    tuple val(label), path("enrichment_${label}"), emit: results
    path "versions.yml", emit: versions, topic: versions

    script:
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
    mkdir -p enrichment_${label} && touch enrichment_${label}/enrichment_summary.csv
    emit_versions.sh "${task.process}" Rscript
    """
}

/*
 * Differential isoform usage. Structure comes from the long-read collapsed GTF, abundances from
 * the short-read quantification — so the GTF must be restricted to the long-read models that
 * survive into THIS reference set, or importRdata's Jaccard check fails.
 */
process ISOFORM_SWITCHING {
    label 'process_high'
    label 'error_ignore'
    tag "$label"

    input:
    path quants, stageAs: 'quants/*'
    tuple val(label), path(fasta)
    path gtf
    path samplesheet

    output:
    tuple val(label), path("isoforms_${label}"), emit: results, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    def rlibs = params.r_libs ? "export R_LIBS_USER=${params.r_libs}" : ''
    """
    ${rlibs}
    grep '^>' ${fasta} | sed 's/^>//;s/[ |].*//' | sort -u > set_ids.txt
    awk -v ids=set_ids.txt '
      BEGIN{ while ((getline l < ids) > 0) keep[l]=1 }
      match(\$0, /transcript_id "[^"]+"/) {
        t = substr(\$0, RSTART+15, RLENGTH-16); if (t in keep) print
      }' ${gtf} > subset.gtf
    echo "[iso] GTF restricted to ${label}: \$(grep -cP '\\ttranscript\\t' subset.gtf) of \$(grep -cP '\\ttranscript\\t' ${gtf}) models"

    mkdir -p isoforms_${label}
    isoform_switching.R --quant-dir quants --fasta ${fasta} --gtf subset.gtf \\
        --samplesheet ${samplesheet} --condition-column '${params.condition_column}' \\
        --outdir isoforms_${label} 2>&1 | tee isoforms_${label}/isoform_switching.log
    emit_versions.sh "${task.process}" Rscript
    """

    stub:
    """
    mkdir -p isoforms_${label} && touch isoforms_${label}/switch_summary.csv
    emit_versions.sh "${task.process}" Rscript
    """
}

/*
 * Gene identity from GENOMIC LOCI rather than from string-stripped transcript ids.
 *
 * String-stripping applies a different rule per assembly source and nothing reconciles gene
 * identity ACROSS sources, so one biological locus recovered by several sources is counted
 * several times. Placing every transcript in one coordinate system and grouping by shared
 * splice junctions (or shared exons) is the stronger lever — at the cost of fusing tandem
 * duplicates and nested genes on a fragmented draft. The report quantifies both directions;
 * read it before adopting the map.
 */
process LOCUS_TX2GENE {
    label 'process_high'
    tag "${label} ${params.locus_groupby} mincov=${params.locus_mincov}"

    input:
    tuple val(label), path(fasta)
    path genome
    tuple val(label2), path(tx2gene)
    tuple val(label3), path(annotation)

    output:
    tuple val("${label}_locus"), path("tx2gene.${label}_locus.tsv")  , emit: tx2gene
    tuple val("${label}_locus"), path("annotation.${label}_locus.tsv"), emit: annotation
    path 'locus_rebuild_report.txt'                                  , emit: report
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    minimap2 -x splice -c -t ${task.cpus} --secondary=no ${genome} ${fasta} > placed.paf 2> minimap2.log
    echo "PAF records: \$(wc -l < placed.paf) | transcripts with an alignment: \$(cut -f1 placed.paf | sort -u | wc -l)"

    # best alignment per transcript -> EXON blocks (+ splice junctions).
    # Merging genomic SPANS would fuse any gene nested in another gene's intron.
    paf_to_exons.py placed.paf ${params.locus_mincov} junctions.tsv > best.bed
    sort -k1,1 -k2,2n -o best.bed best.bed
    echo "placed at >= ${params.locus_mincov} coverage: \$(cut -f4 best.bed | sort -u | wc -l) transcripts, \$(wc -l < best.bed) exon blocks"

    if [ "${params.locus_groupby}" = exon ]; then
        bedtools merge -s -c 4 -o collapse -i best.bed > clusters.bed
        awk -F'\\t' '{n=split(\$NF,a,","); for(i=2;i<=n;i++) print a[1]"\\t"a[i]}' clusters.bed > pairs.tsv
    else
        # multi-exon: link transcripts sharing an identical (scaffold, strand, donor, acceptor)
        sort -k2,2 -k3,3 -k4,4n -k5,5n junctions.tsv \\
          | awk -F'\\t' '{k=\$2"|"\$3"|"\$4"|"\$5; if(k==pk) print first"\\t"\$1; else {pk=k; first=\$1}}' > pairs.tsv
        cut -f1 junctions.tsv | sort -u > multi.ids
        awk -F'\\t' 'NR==FNR{m[\$1];next} !(\$4 in m)' multi.ids best.bed | sort -k1,1 -k2,2n > single.bed
        awk -F'\\t' 'NR==FNR{m[\$1];next}  (\$4 in m)' multi.ids best.bed | sort -k1,1 -k2,2n > multi.bed
        # single-exon transcripts have no junctions: attach them to the locus they sit inside
        if [ -s single.bed ] && [ -s multi.bed ]; then
            bedtools intersect -s -wo -a single.bed -b multi.bed \\
              | awk -F'\\t' -v f=${params.locus_se_overlap} '
                  { ov[\$4"\\t"\$10] += \$13; len[\$4] = \$3-\$2 }
                  END { for (k in ov) { split(k,a,"\\t"); if (ov[k] >= f*len[a[1]]) print a[1]"\\t"a[2] } }' >> pairs.tsv
        fi
        # ...and cluster the leftovers among themselves by overlap
        if [ -s single.bed ]; then
            bedtools merge -s -c 4 -o collapse -i single.bed \\
              | awk -F'\\t' '{n=split(\$NF,a,","); for(i=2;i<=n;i++) print a[1]"\\t"a[i]}' >> pairs.tsv
        fi
    fi
    echo "linking pairs: \$(wc -l < pairs.tsv)"

    locus_tx2gene.py --fasta ${fasta} --pairs pairs.tsv --bed best.bed \\
        --old-tx2gene ${tx2gene} --annotation ${annotation} \\
        --mincov ${params.locus_mincov} --groupby ${params.locus_groupby} \\
        --out-tx2gene tx2gene.${label}_locus.tsv \\
        --out-annotation annotation.${label}_locus.tsv \\
        --report locus_rebuild_report.txt
    emit_versions.sh "${task.process}" minimap2 python3
    """

    stub:
    """
    printf "t1\tLOCUS_000001\n" > tx2gene.${label}_locus.tsv
    touch annotation.${label}_locus.tsv locus_rebuild_report.txt
    emit_versions.sh "${task.process}" minimap2 python3
    """
}
