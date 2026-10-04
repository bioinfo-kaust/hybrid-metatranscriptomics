/*
 * Stage 3 — build the reference transcriptome from four complementary sources:
 *   A  long reads on the genome   (pbmm2 + isoseq collapse)   -> real loci, isoform structure
 *   B  long reads, reference-free (cd-hit-est)                -> recovers what the draft genome misses
 *   C  short reads on the genome  (STAR + StringTie)          -> loci the long reads did not sample
 *   D  short reads de novo        (Trinity)                   -> loci absent from the genome
 * Each source is layered onto the previous one and its contribution is measured, so the ladder
 * is auditable rather than a single opaque merge.
 */

// ---------------------------------------------------------------- Route A: genome-guided long reads
process PBMM2_ALIGN {
    label 'process_high'

    input:
    path seqs
    path genome

    output:
    path 'longreads_aligned.bam'         , emit: bam
    path 'longreads_aligned.flagstat.txt', emit: flagstat
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '--preset ISOSEQ --sort'
    """
    pbmm2 align ${args} -j ${task.cpus} ${genome} ${seqs} longreads_aligned.bam
    samtools flagstat longreads_aligned.bam > longreads_aligned.flagstat.txt
    emit_versions.sh "${task.process}" pbmm2 samtools
    """

    stub:
    """
    touch longreads_aligned.bam longreads_aligned.flagstat.txt
    emit_versions.sh "${task.process}" pbmm2 samtools
    """
}

process ISOSEQ_COLLAPSE {
    label 'process_high'

    input:
    path bam

    output:
    path 'collapsed.gff'        , emit: gff
    path 'collapsed.gtf'        , emit: gtf
    path 'collapsed.fasta'      , emit: fasta
    path 'collapsed.stats.tsv'  , emit: stats
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    isoseq collapse ${args} ${bam} collapsed.gff
    # IsoformSwitchAnalyzeR wants a .gtf; isoseq already writes GFF2-style attributes
    cp collapsed.gff collapsed.gtf
    # `isoseq collapse` decorates headers (PB.1.1|PB.1.1:5118-6219(-)|transcript/569365).
    # Reduce to the bare accession: the gene id is derived by stripping the trailing isoform
    # field, which only works on the bare form — leaving the decoration makes every isoform
    # its own gene and silently turns gene-level DGE into isoform-level.
    awk '/^>/{split(substr(\$0,2),a,"|"); print ">"a[1]; next}{print}' collapsed.fasta > tmp.fa
    mv tmp.fa collapsed.fasta
    seqkit stats -T collapsed.fasta > collapsed.stats.tsv
    cat collapsed.stats.tsv
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch collapsed.gff collapsed.gtf collapsed.fasta collapsed.stats.tsv
    emit_versions.sh "${task.process}" seqkit
    """
}

/*
 * Route A models with GENOME bases: gffread spliced at the collapsed models' exon coordinates.
 * The Iso-Seq consensus carries homopolymer indels that shift the reading frame; the genome at
 * the same coordinates does not (this is what SQANTI3 writes as *_corrected.fasta). Ids and model
 * structure are unchanged, so every downstream map stays valid. Fails if any model is lost.
 */
process LONGREAD_GENOME_SEQ {
    label 'process_low'

    input:
    path gtf
    path genome

    output:
    path 'collapsed.genome.fasta'     , emit: fasta
    path 'collapsed.genome.stats.tsv' , emit: stats
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    gffread -w raw.fa -g ${genome} ${gtf}
    # bare ids (gffread may append attributes) and upper case (soft-masked genomes give lower case)
    seqkit seq -i -u raw.fa > collapsed.genome.fasta
    n_models=\$(awk '\$3=="transcript"' ${gtf} | wc -l)
    n_seq=\$(grep -c '^>' collapsed.genome.fasta)
    [ "\$n_models" = "\$n_seq" ] || { echo "[ERROR] \$n_models models, \$n_seq genome sequences" >&2; exit 1; }
    seqkit stats -T collapsed.genome.fasta > collapsed.genome.stats.tsv
    cat collapsed.genome.stats.tsv
    emit_versions.sh "${task.process}" gffread seqkit
    """

    stub:
    """
    touch collapsed.genome.fasta collapsed.genome.stats.tsv
    emit_versions.sh "${task.process}" gffread seqkit
    """
}

// ---------------------------------------------------------------- clustering helpers
/*
 * Redundancy removal. vsearch rather than cd-hit-est: the identity definition is stated rather
 * than implied (cd-hit measures matches over the SHORTER sequence, which is not obvious from the
 * command line), and --query_cov puts a floor on how much of a sequence must align before it is
 * folded into a representative — cd-hit's -aS equivalent, which this pipeline never set, so a
 * short fragment could be absorbed into a transcript many times its length.
 * --cluster_fast sorts by decreasing length, so the representative is the longest member.
 */
process CLUSTER_NR {
    label 'process_high'
    tag "${label} id=${identity}"

    input:
    tuple val(label), path(fasta)
    val identity

    output:
    tuple val(label), path("${label}.nr.fasta"), emit: fasta
    path "${label}.cluster.log"                , emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: "--iddef 2 --query_cov ${params.cluster_min_cov} --strand both"
    """
    vsearch --cluster_fast ${fasta} --id ${identity} --threads ${task.cpus} ${args} \\
            --centroids ${label}.nr.fasta --uc ${label}.uc \\
            --log ${label}.cluster.log 2>&1
    echo "[${label}] \$(grep -c '^>' ${fasta}) -> \$(grep -c '^>' ${label}.nr.fasta) at id=${identity}" \\
      | tee -a ${label}.cluster.log
    emit_versions.sh "${task.process}" vsearch
    """

    stub:
    """
    touch ${label}.nr.fasta ${label}.cluster.log
    emit_versions.sh "${task.process}" vsearch
    """
}

/*
 * What this source ADDS to a reference set. Alignment-based rather than word-based: a query is
 * already represented when some alignment covers at least `novelty_min_cov` of it at the given
 * identity, and everything else is this source's contribution. minimap2 is already a dependency,
 * the rule is one readable threshold pair instead of cd-hit-est-2d's word-length heuristics, and
 * the per-query evidence stays in the PAF.
 */
process NOVEL_VS_REF {
    label 'process_high'
    tag "${label} vs ${ref_label}"

    input:
    tuple val(label), path(query)
    tuple val(ref_label), path(reference)
    val identity
    val prefix

    output:
    tuple val(label), path("${label}.novel.fasta"), emit: fasta
    path "${label}.novelty.log"                   , emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '-x asm20 -c --secondary=no'
    """
    minimap2 ${args} -t ${task.cpus} ${reference} ${query} > q_vs_ref.paf 2> minimap2.log

    # PAF: 1 qname 2 qlen 3 qstart 4 qend .. 10 matching bases 11 alignment length
    awk -v cov=${params.novelty_min_cov} -v idt=${identity} \\
        '{ qcov = (\$4 - \$3) / \$2; qid = \$10 / \$11; if (qcov >= cov && qid >= idt) print \$1 }' \\
        q_vs_ref.paf | sort -u > represented.ids
    grep '^>' ${query} | sed 's/^>//; s/[[:space:]].*//' | sort -u > all.ids
    comm -23 all.ids represented.ids > novel.ids
    seqkit grep -f novel.ids ${query} > novel.raw.fasta

    # Sanitise the ID ONLY. Rewriting the whole header welds the assembler's description onto the
    # id (">TR_TRINITY_DN1_c0_g1_i1_len_328_path__0_0-327_"), and the tx2gene isoform rule is
    # anchored at the end of the id, so it then never matches and every isoform becomes its own
    # gene. The description carries no information the pipeline uses.
    awk -v p='${prefix}' '/^>/{id=substr(\$1,2); gsub(/[^A-Za-z0-9._-]/,"_",id); print ">" p id; next}{print}' novel.raw.fasta \\
      | seqkit seq -m ${params.min_transcript_len} - > ${label}.novel.fasta

    {
      echo "[${label}] vs ${ref_label} at id>=${identity}, query coverage>=${params.novelty_min_cov}"
      echo "  query sequences      : \$(wc -l < all.ids)"
      echo "  already represented  : \$(wc -l < represented.ids)"
      echo "  novel (>=${params.min_transcript_len} bp): \$(grep -c '^>' ${label}.novel.fasta)"
    } | tee ${label}.novelty.log
    emit_versions.sh "${task.process}" minimap2 seqkit
    """

    stub:
    """
    touch ${label}.novel.fasta ${label}.novelty.log
    emit_versions.sh "${task.process}" minimap2 seqkit
    """
}

// ---------------------------------------------------------------- Route C: genome-guided short reads
process STAR_ALIGN {
    label 'process_high'
    tag "${meta.id}"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam
    path "${meta.id}.Log.final.out"                                    , emit: log
    path "versions.yml", emit: versions, topic: versions

    script:
    // 2-pass, no GTF in the index: junctions are discovered de novo.
    // --outSAMstrandField intronMotif also suppresses spliced alignments whose junctions are all
    // non-canonical and unannotated — conservative, but it drops a few genuine ones.
    def args = task.ext.args ?: '--twopassMode Basic --outSAMstrandField intronMotif'
    """
    STAR --genomeDir ${index} --readFilesIn ${reads[0]} ${reads[1]} --readFilesCommand zcat \\
         --runThreadN ${task.cpus} --outSAMtype BAM SortedByCoordinate \\
         --outFileNamePrefix ${meta.id}. ${args}
    mv ${meta.id}.Aligned.sortedByCoord.out.bam ${meta.id}.bam
    samtools index -@ ${task.cpus} ${meta.id}.bam
    emit_versions.sh "${task.process}" STAR samtools
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai ${meta.id}.Log.final.out
    emit_versions.sh "${task.process}" STAR samtools
    """
}

process STRINGTIE {
    label 'process_medium'
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "${meta.id}.gtf", emit: gtf
    path "versions.yml", emit: versions, topic: versions

    script:
    def strand = params.strandedness == 'reverse' ? '--rf' : (params.strandedness == 'forward' ? '--fr' : '')
    def args   = task.ext.args ?: ''
    """
    stringtie ${bam} -p ${task.cpus} ${strand} ${args} -o ${meta.id}.gtf -l ${meta.id}
    emit_versions.sh "${task.process}" stringtie
    """

    stub:
    """
    touch ${meta.id}.gtf
    emit_versions.sh "${task.process}" stringtie
    """
}

/*
 * Merge the per-sample short-read assemblies and keep only loci the long-read models missed.
 * gffcompare class codes u/x/i/y = intergenic / antisense / intronic / contained-in-intron;
 * everything else is a variant of a locus the long reads already describe.
 */
process STRINGTIE_MERGE {
    label 'process_high'

    input:
    path gtfs, stageAs: 'gtf/*'
    path reference_gff
    path genome
    path reference_ids

    output:
    path 'stringtie_novel.fasta'  , emit: fasta
    path 'stringtie_merged.gtf'   , emit: merged_gtf
    path 'gffcmp.stats'           , emit: stats
    path 'stringtie_novel_ids.txt', emit: novel_ids
    path "versions.yml", emit: versions, topic: versions

    script:
    def codes = params.stringtie_class_codes
    def has_ref = reference_gff.name != 'NO_FILE_GFF'
    def has_ids = reference_ids.name != 'NO_FILE'
    """
    ls gtf/* > gtf_list.txt
    stringtie --merge -p ${task.cpus} -o stringtie_merged.gtf gtf_list.txt

    if ${has_ref}; then
        # Compare against the long-read models that SURVIVED the contaminant audit. Using the raw
        # collapse output meant a StringTie locus matching a model later deleted as contaminant was
        # judged "not novel" and disappeared from every rung.
        REF=${reference_gff}
        if ${has_ids}; then
            awk 'NR==FNR{keep[\$1]=1; next}
                 { if (match(\$0, /transcript_id "[^"]+"/)) {
                     t = substr(\$0, RSTART+15, RLENGTH-16); if (t in keep) print }
                 }' ${reference_ids} ${reference_gff} > reference.audited.gff
            REF=reference.audited.gff
            echo "reference models: \$(grep -c 'transcript_id' ${reference_gff}) -> \$(grep -c 'transcript_id' reference.audited.gff) after the audit"
        fi

        gffcompare -r \$REF -o gffcmp stringtie_merged.gtf > gffcompare.log 2>&1 || true
        # POSIX awk: match() takes two arguments, and the 3-argument gawk form used to fail
        # silently here (behind || true), leaving an empty StringTie contribution that looked
        # like a real result. RSTART/RLENGTH carry the same information portably.
        awk -v codes="${codes}" '
            \$3 == "transcript" {
                cc = ""; tid = ""
                if (match(\$0, /class_code "[^"]+"/)) cc  = substr(\$0, RSTART+12, RLENGTH-13)
                if (match(\$0, /transcript_id "[^"]+"/)) tid = substr(\$0, RSTART+15, RLENGTH-16)
                if (cc != "" && tid != "" && cc ~ "^(" codes ")\$") print tid
            }' gffcmp.annotated.gtf | sort -u > stringtie_novel_ids.txt
        # select records by PARSED transcript_id: -Fw matched the id as a substring anywhere on
        # the line, so a gene_id sharing the prefix pulled in records that were never selected.
        awk 'NR==FNR{keep[\$1]=1; next}
             { if (match(\$0, /transcript_id "[^"]+"/)) {
                 t = substr(\$0, RSTART+15, RLENGTH-16); if (t in keep) print }
             }' stringtie_novel_ids.txt gffcmp.annotated.gtf > novel.gtf
    else
        # no long-read models to compare against: every merged transcript is the base set
        touch gffcmp.stats
        grep -oE 'transcript_id "[^"]+"' stringtie_merged.gtf | sed 's/transcript_id "//;s/"//' \\
          | sort -u > stringtie_novel_ids.txt
        cp stringtie_merged.gtf novel.gtf
    fi
    echo "short-read transcripts taken forward: \$(wc -l < stringtie_novel_ids.txt)"
    gffread -w novel.raw.fasta -g ${genome} novel.gtf
    awk '/^>/{sub(/^>/,">${params.prefix_stringtie}")}1' novel.raw.fasta \\
      | seqkit seq -m ${params.min_transcript_len} - > stringtie_novel.fasta
    echo "novel >=${params.min_transcript_len} bp: \$(grep -c '^>' stringtie_novel.fasta)"
    emit_versions.sh "${task.process}" stringtie gffcompare gffread
    """

    stub:
    """
    touch stringtie_novel.fasta stringtie_merged.gtf gffcmp.stats stringtie_novel_ids.txt
    emit_versions.sh "${task.process}" stringtie gffcompare gffread
    """
}

// ---------------------------------------------------------------- Route D: de novo short reads
process TRINITY {
    label 'process_high'
    label 'process_long'
    label 'process_high_memory'

    input:
    path reads1, stageAs: 'r1/*'
    path reads2, stageAs: 'r2/*'

    output:
    path 'trinity.Trinity.fasta', emit: fasta
    path 'trinity.stats.tsv'    , emit: stats
    path "versions.yml", emit: versions, topic: versions

    script:
    def ss   = params.strandedness == 'reverse' ? '--SS_lib_type RF' : (params.strandedness == 'forward' ? '--SS_lib_type FR' : '')
    def args = task.ext.args ?: '--no_version_check'
    """
    L=\$(ls r1/* | paste -sd, -)
    R=\$(ls r2/* | paste -sd, -)
    Trinity --seqType fq --max_memory ${params.trinity_max_memory} --CPU ${task.cpus} \\
            ${ss} --left "\$L" --right "\$R" --output trinity ${args}
    # Trinity writes <outdir>.Trinity.fasta next to the output dir, but the exact location has
    # moved between releases — resolve it rather than assuming
    if [ ! -s trinity.Trinity.fasta ]; then
        SRC=\$(ls trinity/Trinity.fasta trinity.Trinity.fasta trinity/trinity.Trinity.fasta 2>/dev/null | head -1)
        [ -n "\$SRC" ] || { echo "[ERROR] Trinity produced no assembly" >&2; exit 1; }
        cp "\$SRC" trinity.Trinity.fasta
    fi
    seqkit stats -T trinity.Trinity.fasta > trinity.stats.tsv
    cat trinity.stats.tsv
    emit_versions.sh "${task.process}" Trinity seqkit
    """

    stub:
    """
    touch trinity.Trinity.fasta trinity.stats.tsv
    emit_versions.sh "${task.process}" Trinity seqkit
    """
}

// ---------------------------------------------------------------- ladder assembly
/* Concatenate a base set with the novel contributions of one or more sources. */
process MERGE_SET {
    label 'process_low'
    tag "$label"

    input:
    tuple val(label), path(parts, stageAs: 'part/*')

    output:
    tuple val(label), path("${label}.fasta"), emit: fasta
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    cat part/* > ${label}.fasta
    echo "[${label}] transcripts: \$(grep -c '^>' ${label}.fasta)"
    emit_versions.sh "${task.process}" seqkit
    """

    stub:
    """
    touch ${label}.fasta
    emit_versions.sh "${task.process}" seqkit
    """
}

process MAKE_TX2GENE {
    label 'process_single'
    tag "$label"

    input:
    tuple val(label), path(fasta)

    output:
    tuple val(label), path("tx2gene.${label}.tsv"), emit: tx2gene
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    make_tx2gene.py --fasta ${fasta} --rules '${params.tx2gene_rules}' --out tx2gene.${label}.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    printf "t1\tg1\n" > tx2gene.${label}.tsv
    emit_versions.sh "${task.process}" python3
    """
}

/* Restrict a tx2gene map to the transcripts that survive into a deduplicated set. */
process SUBSET_TX2GENE {
    label 'process_single'
    tag "$label"

    input:
    tuple val(label), path(fasta)
    path tx2gene

    output:
    tuple val(label), path("tx2gene.${label}.tsv"), emit: tx2gene
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    grep '^>' ${fasta} | sed 's/^>//; s/[ |].*//' | sort -u > ids.txt
    awk 'NR==FNR{k[\$1];next} (\$1 in k)' ids.txt ${tx2gene} > tx2gene.${label}.tsv
    echo "[${label}] transcripts \$(wc -l < ids.txt), genes \$(cut -f2 tx2gene.${label}.tsv | sort -u | wc -l)"
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    printf "t1\tg1\n" > tx2gene.${label}.tsv
    emit_versions.sh "${task.process}" python3
    """
}

/* The auditable ladder: one row per reference variant with transcripts, genes and BUSCO. */
process ROUTE_CONTRIBUTION {
    label 'process_single'

    input:
    path stats  , stageAs: 'stats/*'
    path tx2gene, stageAs: 't2g/*'
    path busco  , stageAs: 'busco/*'

    output:
    path 'route_contribution.tsv', emit: table
    path "versions.yml", emit: versions, topic: versions

    script:
    """
    route_contribution.py --stats stats --tx2gene t2g --busco busco --out route_contribution.tsv
    column -t -s\$'\\t' route_contribution.tsv
    emit_versions.sh "${task.process}" python3
    """

    stub:
    """
    touch route_contribution.tsv
    emit_versions.sh "${task.process}" python3
    """
}

process BUSCO {
    label 'process_high'
    label 'error_ignore'
    tag "$label"

    input:
    tuple val(label), path(fasta)

    output:
    tuple val(label), path("${label}.busco.txt"), emit: summary
    path "busco_${label}/**"                    , emit: dir, optional: true
    path "versions.yml", emit: versions, topic: versions

    script:
    def dl = params.busco_download_path ? "--download_path ${params.busco_download_path}" : ''
    """
    unset PYTHONPATH
    # Serial ids: only the completeness string is consumed downstream, and BUSCO crashes
    # (KeyError in GenomeAnalysis.get_exon_records) when metaeuk truncates very long de novo
    # headers into colliding gene ids.
    awk '/^>/{printf(">t%d\\n", ++i); next}{print}' ${fasta} > in.fa
    busco -i in.fa -l ${params.busco_lineage} -m ${params.busco_mode} -c ${task.cpus} \\
          --out_path . -o busco_${label} ${dl} -f > busco.log 2>&1 || echo "[WARN] busco exited non-zero for ${label}" >&2
    grep -hoE "C:[0-9.]+%\\[S:[0-9.]+%,D:[0-9.]+%\\],F:[0-9.]+%,M:[0-9.]+%" \\
        busco_${label}/short_summary.*.txt 2>/dev/null | head -1 > ${label}.busco.txt || true
    [ -s ${label}.busco.txt ] || echo NA > ${label}.busco.txt
    echo "[${label}] \$(cat ${label}.busco.txt)"
    emit_versions.sh "${task.process}" busco
    """

    stub:
    """
    echo "C:70.0%[S:50.0%,D:20.0%],F:10.0%,M:20.0%" > ${label}.busco.txt
    emit_versions.sh "${task.process}" busco
    """
}

/*
 * SQANTI3 structural QC of the long-read models against the short-read annotation.
 *
 * It was in the project plan from the start as the isoform-quality gate and never ran. The
 * categories (FSM / ISM / NIC / NNC / antisense / intergenic) and the RT-switching and
 * intra-priming flags say which collapse products are supported by an independent assembly and
 * which look like artefacts. This reports; it does not yet filter the reference.
 */
process SQANTI3_QC {
    label 'process_high'
    label 'error_ignore'
    tag "$label"

    input:
    tuple val(label), path(isoforms_gtf)
    path reference_gtf
    path genome

    output:
    path "sqanti3_${label}"            , emit: results, optional: true
    path "sqanti3_${label}.summary.txt", emit: summary
    path "versions.yml", emit: versions, topic: versions

    script:
    def args = task.ext.args ?: '--skipORF --report skip'
    """
    # The module ships its own python env, but ~/.local/lib/python3.9 shadows it and pulls in a
    # pybedtools built against a newer libstdc++ than the module carries. Ignore user site-packages.
    export PYTHONNOUSERSITE=1

    sqanti3_qc.py ${isoforms_gtf} ${reference_gtf} ${genome} \\
        -g -t ${task.cpus} -d sqanti3_${label} -o ${label} ${args} > sqanti3.log 2>&1 || true

    CLS=\$(ls sqanti3_${label}/*classification.txt 2>/dev/null | head -1)
    {
      echo "SQANTI3 structural categories — ${label} vs the short-read annotation"
      if [ -n "\$CLS" ] && [ -s "\$CLS" ]; then
          awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) h[\$i]=i; next}
                       { c[\$h["structural_category"]]++; n++ }
                       END { for (k in c) printf "  %-28s %8d  %5.1f%%\\n", k, c[k], 100*c[k]/n;
                             printf "  %-28s %8d\\n", "total", n }' "\$CLS" | sort -k2,2nr
      else
          echo "  [WARN] SQANTI3 produced no classification table — see sqanti3.log"
      fi
    } | tee sqanti3_${label}.summary.txt
    emit_versions.sh "${task.process}" sqanti3_qc.py
    """

    stub:
    """
    mkdir -p sqanti3_${label} && touch sqanti3_${label}.summary.txt
    emit_versions.sh "${task.process}" sqanti3_qc.py
    """
}
