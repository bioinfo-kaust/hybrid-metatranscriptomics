/*
 * Stage 3 — build the reference transcriptome as an auditable ladder.
 *
 * Each source is layered onto the previous one and only its NOVEL contribution is added, so the
 * table at the end says what every source was worth. The rungs:
 *
 *   isoseq                           Iso-Seq models on the genome (Route A) — real loci, isoforms
 *   isoseq_stringtie                 + genome-guided short-read loci the long reads missed
 *   isoseq_stringtie_trinity         + de novo short-read loci              (Trinity)
 *   isoseq_stringtie_lrfree          + reference-free long-read loci only   (vsearch)
 *   isoseq_stringtie_trinity_lrfree  all four sources merged
 *   all_sources_dedup                cross-source dedup  <- the assembly: annotated, BUSCO-scored,
 *                                                          then expression-filtered for quantification
 *
 * Each name reads as the sources it contains, in the order they were layered, so a rung says
 * what is in it rather than when it was made. The expression filter derives one more rung,
 * all_sources_dedup_expressed, which is the set actually quantified and tested.
 *
 * Without a genome there is no Route A or StringTie, so the isoseq_stringtie rung is empty: the
 * de novo and reference-free sets are then novel against nothing, and the later rungs hold only
 * them. The names stay the same so dge_reference and the reports need no special case.
 *
 * The lean set exists because the other rungs count one biological locus several times when
 * several sources recover it. Deduplication reduces that; it does not eliminate it, which is
 * what the locus rebuild in stage 6 addresses.
 */

include { PBMM2_ALIGN                                } from '../../../modules/nf-core/pbmm2/align'
include { SAMTOOLS_FLAGSTAT                          } from '../../../modules/nf-core/samtools/flagstat'
include { ISOSEQ_COLLAPSE                            } from '../../../modules/local/isoseq_collapse'
include { LONGREAD_GENOME_SEQ                        } from '../../../modules/local/longread_genome_seq'
include { STAR_ALIGN                                 } from '../../../modules/nf-core/star/align'
include { STRINGTIE_STRINGTIE                        } from '../../../modules/nf-core/stringtie/stringtie'
include { STRINGTIE_MERGE                            } from '../../../modules/nf-core/stringtie/merge'
include { FILTER_GTF_BY_IDS                          } from '../../../modules/local/filter_gtf_by_ids'
include { GFFCOMPARE                                 } from '../../../modules/nf-core/gffcompare'
include { STRINGTIE_NOVEL                            } from '../../../modules/local/stringtie_novel'
include { TRINITY                                    } from '../../../modules/nf-core/trinity'
include { MERGE_SET as MERGE_ISOSEQ                  } from '../../../modules/local/merge_set'
include { MERGE_SET as MERGE_DEDUP                   } from '../../../modules/local/merge_set'
include { MERGE_SET as MERGE_ISOSEQ_STRINGTIE        } from '../../../modules/local/merge_set'
include { MERGE_SET as MERGE_WITH_TRINITY            } from '../../../modules/local/merge_set'
include { MERGE_SET as MERGE_WITH_LRFREE             } from '../../../modules/local/merge_set'
include { MERGE_SET as MERGE_ALL_SOURCES             } from '../../../modules/local/merge_set'
include { MAKE_TX2GENE                               } from '../../../modules/local/make_tx2gene'
include { SUBSET_TX2GENE                             } from '../../../modules/local/subset_tx2gene'
include { ROUTE_CONTRIBUTION                         } from '../../../modules/local/route_contribution'
include { SEQKIT_STATS as SEQ_STATS                  } from '../../../modules/nf-core/seqkit/stats'
include { CLASSIFY_CONTAM as AUDIT_LONGREAD          } from '../../../modules/local/classify_contam'
include { CLASSIFY_CONTAM as AUDIT_LRFREE            } from '../../../modules/local/classify_contam'
include { CLASSIFY_CONTAM as AUDIT_DENOVO            } from '../../../modules/local/classify_contam'
include { VSEARCH_CLUSTER as CLUSTER_LRFREE          } from '../../../modules/nf-core/vsearch/cluster'
include { VSEARCH_CLUSTER as CLUSTER_LRFREE_NR       } from '../../../modules/nf-core/vsearch/cluster'
include { VSEARCH_CLUSTER as CLUSTER_DENOVO          } from '../../../modules/nf-core/vsearch/cluster'
include { VSEARCH_CLUSTER as CLUSTER_LEAN            } from '../../../modules/nf-core/vsearch/cluster'
include { NOVEL_VS_REF as NOVEL_DENOVO               } from '../../../modules/local/novel_vs_ref'
include { NOVEL_VS_REF as NOVEL_LRFREE_A             } from '../../../modules/local/novel_vs_ref'
include { NOVEL_VS_REF as NOVEL_LRFREE_B             } from '../../../modules/local/novel_vs_ref'
include { SQANTI3_QC                                 } from '../../../modules/local/sqanti3_qc'
include { TRANSCRIPTOME_BUSCO                        } from '../transcriptome_busco'
include { runLongReads; runGenomeGuided; runStringtie; runSqanti; runDecontam } from '../utils_nfcore_hybridmetatranscriptomics_pipeline'

workflow BUILD_TRANSCRIPTOME {

    take:
    ch_longreads      // [ meta, fasta ] or empty
    ch_reads          // [ meta, [r1, r2] ]  the experiment's libraries plus any assembly-only ones
    ch_genome         // path (value)
    ch_contaminant    // path (value)
    ch_star_index     // [ meta, dir ] (value)

    main:
    ch_sets   = channel.empty()      // [ meta, fasta ] for every ladder rung
    ch_counts = channel.empty()      // audit reports
    ch_qc     = channel.empty()

    def do_long   = runLongReads()
    def do_routea = do_long && runGenomeGuided()
    def do_routeb = do_long && !params.skip_reffree
    def do_st     = runStringtie()
    def has_base  = do_routea || do_st     // anything for the later routes to be novel against
    def audit     = runDecontam()

    // The residual audit drops ONLY what positively places on a contaminant. Read-level removal
    // already keeps reads that align nowhere — "unmapped is not evidence of contamination" on a
    // fragmented draft genome — and the audit contradicted that by deleting every assembled
    // sequence that aligned to neither reference. That is the de novo route's whole purpose, so
    // it keeps them and reports the count.
    def audit_keep = 'host,ambiguous,unmapped'

    ch_genome_meta = ch_genome.map { fa -> [ [id: 'genome'], fa ] }

    // ---------------------------------------------------- Route A: long reads on the genome
    ch_collapsed_gff = channel.empty()
    ch_collapsed_gtf = channel.empty()
    ch_base          = channel.empty()

    if (do_routea) {
        PBMM2_ALIGN(ch_longreads, ch_genome_meta)
        SAMTOOLS_FLAGSTAT(PBMM2_ALIGN.out.bam.join(PBMM2_ALIGN.out.index, remainder: true)
            .map { meta, bam, bai -> [ meta, bam, bai ?: [] ] })
        ISOSEQ_COLLAPSE(PBMM2_ALIGN.out.bam)
        ch_collapsed_gff = ISOSEQ_COLLAPSE.out.gff
        ch_collapsed_gtf = ISOSEQ_COLLAPSE.out.gtf

        // genome bases at the model coordinates, not the consensus: see LONGREAD_GENOME_SEQ
        ch_routea_fa = ISOSEQ_COLLAPSE.out.fasta
        if (params.longread_genome_sequence) {
            LONGREAD_GENOME_SEQ(ISOSEQ_COLLAPSE.out.gtf, ch_genome)
            ch_routea_fa = LONGREAD_GENOME_SEQ.out.fasta
        }
        ch_routea = ch_routea_fa.map { _meta, fa -> [ [id: 'routeA'], fa ] }
        if (audit) {
            // residual audit: read-level removal already happened, so this is the QC gate
            AUDIT_LONGREAD(ch_routea, ch_genome, ch_contaminant, audit_keep)
            ch_counts = ch_counts.mix(AUDIT_LONGREAD.out.counts)
            ch_routea = AUDIT_LONGREAD.out.fasta
            // gffcompare is scored against the models that SURVIVED the audit: against the raw
            // collapse output, a StringTie locus matching a model later deleted as contaminant was
            // judged "not novel" and disappeared from every rung
            ch_collapsed_gff = FILTER_GTF_BY_IDS(ISOSEQ_COLLAPSE.out.gff, AUDIT_LONGREAD.out.host_ids).gtf
        }
        // every ladder rung is materialised as <label>.fasta so the reference set is a file,
        // not something the reader has to reconstruct from a work directory
        ch_base = MERGE_ISOSEQ(ch_routea.map { _meta, fa -> [ [id: 'isoseq'], [fa] ] }).fasta
        ch_sets = ch_sets.mix(ch_base)
    }

    // ------------------------------------------ Route C: genome-guided short-read assembly
    ch_isoseq_stringtie = ch_base
    ch_merged_gtf       = channel.empty()
    if (do_st) {
        STAR_ALIGN(ch_reads, ch_star_index, [ [:], [] ], true)
        ch_qc = ch_qc.mix(STAR_ALIGN.out.log_final.map { _meta, log -> log })
        STRINGTIE_STRINGTIE(STAR_ALIGN.out.bam_sorted_aligned.map { meta, bam -> [ meta, bam, [] ] }, '', [])

        // merged in file-name order, so the merge does not depend on task completion order
        STRINGTIE_MERGE(
            STRINGTIE_STRINGTIE.out.transcript_gtf
                .map { _meta, gtf -> gtf }
                .toSortedList { a, b -> a.name <=> b.name }
                .map { gtfs -> [ [id: 'stringtie_merged'], gtfs ] },
            [ [:], [] ])
        ch_merged_gtf = STRINGTIE_MERGE.out.merged_gtf

        // with long-read models, keep only the loci they missed; without them every merged
        // transcript is taken
        if (do_routea) {
            GFFCOMPARE(ch_merged_gtf, [ [:], [], [] ], ch_collapsed_gff)
            ch_annotated = GFFCOMPARE.out.annotated_gtf.map { _meta, gtf -> gtf }
        } else {
            ch_annotated = channel.value(file("${projectDir}/assets/NO_FILE_GFFCMP", checkIfExists: true))
        }
        STRINGTIE_NOVEL(ch_merged_gtf, ch_annotated, ch_genome)
        ch_st_novel = STRINGTIE_NOVEL.out.fasta.map { _meta, fa -> fa }

        // with long reads present the short-read set only contributes the loci they missed;
        // without them it IS the base set
        ch_boosted_parts = do_routea
            ? ch_base.map { _meta, fa -> fa }.combine(ch_st_novel).map { base, novel -> [ [id: 'isoseq_stringtie'], [base, novel] ] }
            : ch_st_novel.map { novel -> [ [id: 'isoseq_stringtie'], [novel] ] }
        ch_isoseq_stringtie = MERGE_ISOSEQ_STRINGTIE(ch_boosted_parts).fasta
        ch_sets = ch_sets.mix(ch_isoseq_stringtie)
    }

    // Structural QC of the long-read models against an independent short-read annotation.
    // Reports only — the categories say which collapse products a second assembly supports.
    ch_sqanti = channel.empty()
    if (do_routea && do_st && runSqanti()) {
        SQANTI3_QC(ch_collapsed_gtf.map { _meta, gtf -> [ [id: 'isoseq'], gtf ] }, ch_merged_gtf, ch_genome)
        ch_sqanti = SQANTI3_QC.out.summary
    }

    // what the de novo and reference-free sets are compared against
    ch_novelty_ref = has_base ? ch_isoseq_stringtie.first()
                              : channel.value([ [id: 'none'], file("${projectDir}/assets/NO_FILE_REF", checkIfExists: true) ])

    // -------------------------------------------------- Route D: de novo short-read assembly
    ch_denovo_novel = channel.empty()
    if (!params.skip_denovo) {
        // every library in R1 file-name order, interleaved R1/R2 per library
        TRINITY(ch_reads
            .toSortedList { a, b -> a[1][0].name <=> b[1][0].name }
            .map { libs -> [ [id: 'trinity'], libs.collect { _meta, reads -> reads }.flatten() ] })

        ch_trin = TRINITY.out.transcript_fasta.map { _meta, fa -> [ [id: 'denovo'], fa ] }
        if (audit) {
            AUDIT_DENOVO(ch_trin, ch_genome, ch_contaminant, audit_keep)
            ch_trin = AUDIT_DENOVO.out.fasta
            ch_counts = ch_counts.mix(AUDIT_DENOVO.out.counts)
        }
        // collapse the assembler's own isoform redundancy, then keep only what the base lacks
        CLUSTER_DENOVO(ch_trin)
        NOVEL_DENOVO(CLUSTER_DENOVO.out.centroids, ch_novelty_ref, params.cdhit_dedup_id, params.prefix_denovo)
        ch_denovo_novel = NOVEL_DENOVO.out.fasta.map { _meta, fa -> fa }
    }

    // ------------------------------------------- Route B: reference-free long-read collapse
    ch_lrfree_novel_vs_base = channel.empty()
    ch_lrfree_nr            = channel.empty()
    if (do_routeb) {
        CLUSTER_LRFREE(ch_longreads.map { _meta, fa -> [ [id: 'lrfree'], fa ] })
        ch_lrfree = CLUSTER_LRFREE.out.centroids
        if (audit) {
            AUDIT_LRFREE(ch_lrfree, ch_genome, ch_contaminant, audit_keep)
            ch_lrfree = AUDIT_LRFREE.out.fasta
            ch_counts = ch_counts.mix(AUDIT_LRFREE.out.counts)
        }
        CLUSTER_LRFREE_NR(ch_lrfree)
        ch_lrfree_nr = CLUSTER_LRFREE_NR.out.centroids

        // (a) what Route B adds to the base, on its own
        NOVEL_LRFREE_A(ch_lrfree_nr, ch_novelty_ref, params.cdhit_dedup_id, params.prefix_lrfree)
        ch_lrfree_novel_vs_base = NOVEL_LRFREE_A.out.fasta.map { _meta, fa -> fa }
    }

    // ---------------------------------------------------------------- ladder assembly
    ch_with_trinity = channel.empty()
    if (!params.skip_denovo) {
        ch_with_trinity = MERGE_WITH_TRINITY(has_base
            ? ch_isoseq_stringtie.map { _meta, fa -> fa }.combine(ch_denovo_novel)
                .map { base, novel -> [ [id: 'isoseq_stringtie_trinity'], [base, novel] ] }
            : ch_denovo_novel.map { novel -> [ [id: 'isoseq_stringtie_trinity'], [novel] ] }).fasta
        ch_sets = ch_sets.mix(ch_with_trinity)
    }

    if (do_routeb) {
        ch_with_lrfree = MERGE_WITH_LRFREE(has_base
            ? ch_isoseq_stringtie.map { _meta, fa -> fa }.combine(ch_lrfree_novel_vs_base)
                .map { base, novel -> [ [id: 'isoseq_stringtie_lrfree'], [base, novel] ] }
            : ch_lrfree_novel_vs_base.map { novel -> [ [id: 'isoseq_stringtie_lrfree'], [novel] ] }).fasta
        ch_sets = ch_sets.mix(ch_with_lrfree)
    }

    // (b) what Route B still adds once the de novo set is already in — its unique contribution
    if (do_routeb && !params.skip_denovo) {
        NOVEL_LRFREE_B(ch_lrfree_nr, ch_with_trinity, params.cdhit_dedup_id, params.prefix_lrfree)
        ch_all_sources = MERGE_ALL_SOURCES(
            ch_with_trinity.map { _meta, fa -> fa }
                .combine(NOVEL_LRFREE_B.out.fasta.map { _meta, fa -> fa })
                .map { base, novel -> [ [id: 'isoseq_stringtie_trinity_lrfree'], [base, novel] ] }).fasta
        ch_sets = ch_sets.mix(ch_all_sources)
    }

    // the widest set available becomes the input to the cross-source deduplication
    def widest_label = params.skip_denovo ? (do_routeb ? 'isoseq_stringtie_lrfree' : 'isoseq_stringtie')
                                          : (do_routeb ? 'isoseq_stringtie_trinity_lrfree' : 'isoseq_stringtie_trinity')
    ch_widest = ch_sets.filter { meta, _fa -> meta.id == widest_label }

    CLUSTER_LEAN(ch_widest.map { _meta, fa -> [ [id: 'dedup'], fa ] })
    ch_lean = MERGE_DEDUP(CLUSTER_LEAN.out.centroids.map { _meta, fa -> [ [id: 'all_sources_dedup'], [fa] ] }).fasta

    // ---------------------------------------------------------------- maps, stats, BUSCO
    // The lean set's map is the widest set's map restricted to the surviving transcripts, so
    // gene identity stays consistent between the two.
    MAKE_TX2GENE(ch_sets)
    SUBSET_TX2GENE(ch_lean, MAKE_TX2GENE.out.tx2gene.filter { meta, _t -> meta.id == widest_label }.map { _meta, t -> t })

    ch_all_sets = ch_sets.mix(ch_lean)
    ch_tx2gene  = MAKE_TX2GENE.out.tx2gene.mix(SUBSET_TX2GENE.out.tx2gene)

    SEQ_STATS(ch_all_sets)
    ch_busco = channel.empty()
    if (!params.skip_busco) {
        TRANSCRIPTOME_BUSCO(ch_all_sets)
        ch_busco = TRANSCRIPTOME_BUSCO.out.summary
        ch_qc    = ch_qc.mix(TRANSCRIPTOME_BUSCO.out.mqc)
    }

    ROUTE_CONTRIBUTION(
        SEQ_STATS.out.stats.map { _meta, f -> f }.collect(),
        ch_tx2gene.map { _meta, f -> f }.collect(),
        params.skip_busco ? channel.fromPath("${projectDir}/assets/NO_FILE_BUSCO").collect()
                          : ch_busco.map { _meta, f -> f }.collect())

    // ------------------------------------------------------- pick the quantification reference
    ch_selected     = ch_all_sets.filter { meta, _fa -> meta.id == params.dge_reference }
    ch_selected_t2g = ch_tx2gene.filter { meta, _t -> meta.id == params.dge_reference }

    emit:
    sets             = ch_all_sets       // [ meta, fasta ]
    tx2gene          = ch_tx2gene        // [ meta, tsv ]
    selected         = ch_selected       // [ meta, fasta ]
    selected_tx2gene = ch_selected_t2g   // [ meta, tsv ]
    collapsed_gtf    = ch_collapsed_gtf  // [ meta, gtf ]
    ladder           = ROUTE_CONTRIBUTION.out.table
    busco            = ch_busco          // [ meta, <id>.busco.txt ]
    sqanti           = ch_sqanti
    counts           = ch_counts
    qc               = ch_qc
}
