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
    tag "${meta.id} ${params.locus_groupby} mincov=${params.locus_mincov}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/minimap2_bedtools_python_gawk:LOCUS_SIF_PENDING'
        : 'community.wave.seqera.io/library/minimap2_bedtools_python_gawk:LOCUS_DOCKER_PENDING' }"

    input:
    tuple val(meta), path(fasta)
    path genome
    tuple val(meta2), path(tx2gene)
    tuple val(meta3), path(annotation)

    output:
    tuple val(meta_out), path("tx2gene.${meta.id}_locus.tsv")   , emit: tx2gene
    tuple val(meta_out), path("annotation.${meta.id}_locus.tsv"), emit: annotation
    path 'locus_rebuild_report.txt'                             , emit: report
    path "versions.yml"                                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    meta_out  = meta + [id: "${meta.id}_locus"]
    def label = meta.id
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
    emit_versions.sh "${task.process}" minimap2 bedtools python3
    """

    stub:
    meta_out  = meta + [id: "${meta.id}_locus"]
    """
    printf "t1\\tLOCUS_000001\\n" > tx2gene.${meta.id}_locus.tsv
    touch annotation.${meta.id}_locus.tsv locus_rebuild_report.txt
    emit_versions.sh "${task.process}" minimap2 bedtools python3
    """
}
