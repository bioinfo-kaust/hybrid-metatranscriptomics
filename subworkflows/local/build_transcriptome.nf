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
 * The lean set exists because the other rungs count one biological locus several times when
 * several sources recover it. Deduplication reduces that; it does not eliminate it, which is
 * what the locus rebuild in stage 6 addresses.
 */

include { PBMM2_ALIGN     } from '../../modules/local/assembly'
include { ISOSEQ_COLLAPSE } from '../../modules/local/assembly'
include { LONGREAD_GENOME_SEQ } from '../../modules/local/assembly'
include { STAR_ALIGN      } from '../../modules/local/assembly'
include { STRINGTIE       } from '../../modules/local/assembly'
include { STRINGTIE_MERGE } from '../../modules/local/assembly'
include { TRINITY         } from '../../modules/local/assembly'
include { MERGE_SET as MERGE_ISOSEQ           } from '../../modules/local/assembly'
include { MERGE_SET as MERGE_DEDUP            } from '../../modules/local/assembly'
include { MERGE_SET as MERGE_ISOSEQ_STRINGTIE } from '../../modules/local/assembly'
include { MERGE_SET as MERGE_WITH_TRINITY     } from '../../modules/local/assembly'
include { MERGE_SET as MERGE_WITH_LRFREE      } from '../../modules/local/assembly'
include { MERGE_SET as MERGE_ALL_SOURCES      } from '../../modules/local/assembly'
include { MAKE_TX2GENE    } from '../../modules/local/assembly'
include { SUBSET_TX2GENE  } from '../../modules/local/assembly'
include { ROUTE_CONTRIBUTION } from '../../modules/local/assembly'
include { BUSCO           } from '../../modules/local/assembly'
include { SEQ_STATS       } from '../../modules/local/references'
include { CLASSIFY_CONTAM as AUDIT_LONGREAD } from '../../modules/local/decontam'
include { CLASSIFY_CONTAM as AUDIT_LRFREE   } from '../../modules/local/decontam'
include { CLASSIFY_CONTAM as AUDIT_DENOVO   } from '../../modules/local/decontam'
include { CLUSTER_NR as CLUSTER_LRFREE    } from '../../modules/local/assembly'
include { CLUSTER_NR as CLUSTER_LRFREE_NR } from '../../modules/local/assembly'
include { CLUSTER_NR as CLUSTER_DENOVO    } from '../../modules/local/assembly'
include { CLUSTER_NR as CLUSTER_LEAN      } from '../../modules/local/assembly'
include { NOVEL_VS_REF as NOVEL_DENOVO   } from '../../modules/local/assembly'
include { NOVEL_VS_REF as NOVEL_LRFREE_A } from '../../modules/local/assembly'
include { NOVEL_VS_REF as NOVEL_LRFREE_B } from '../../modules/local/assembly'
include { SQANTI3_QC      } from '../../modules/local/assembly'

// The residual audit drops ONLY what positively places on a contaminant. Read-level removal
// already keeps reads that align nowhere — "unmapped is not evidence of contamination" on a
// fragmented draft genome — and the audit contradicted that by deleting every assembled sequence
// that aligned to neither reference. That is the de novo route's whole purpose, so it keeps them
// and reports the count.
def AUDIT_KEEP = 'host,ambiguous,unmapped'

workflow BUILD_TRANSCRIPTOME {

    take:
    ch_longreads      // [ 'longreads', fasta ] or empty
    ch_reads          // [ meta, [r1, r2] ]  the experiment's libraries plus any assembly-only ones
    ch_genome
    ch_contaminant
    ch_star_index

    main:
    ch_sets   = Channel.empty()      // [ label, fasta ] for every ladder rung
    ch_counts = Channel.empty()      // audit reports
    ch_qc     = Channel.empty()

    def do_long   = !params.skip_longread
    def do_routea = do_long && !params.skip_genome_guided
    def do_routeb = do_long && !params.skip_reffree
    def audit     = !params.skip_decontam

    // ---------------------------------------------------- Route A: long reads on the genome
    ch_collapsed_gff = Channel.fromPath("${projectDir}/assets/NO_FILE_GFF", checkIfExists: false)
    // ids of the long-read models that survived the contaminant audit; gffcompare is scored
    // against these, not against the raw collapse output
    ch_routea_ids = Channel.fromPath("${projectDir}/assets/NO_FILE").first()
    ch_collapsed_gtf = Channel.empty()
    ch_base = Channel.empty()

    if (do_routea) {
        PBMM2_ALIGN(ch_longreads.map { l, fa -> fa }, ch_genome)
        ISOSEQ_COLLAPSE(PBMM2_ALIGN.out.bam)
        ch_collapsed_gff = ISOSEQ_COLLAPSE.out.gff
        ch_collapsed_gtf = ISOSEQ_COLLAPSE.out.gtf

        // genome bases at the model coordinates, not the consensus: see LONGREAD_GENOME_SEQ
        ch_routea_fa = ISOSEQ_COLLAPSE.out.fasta
        if (params.longread_genome_sequence) {
            LONGREAD_GENOME_SEQ(ISOSEQ_COLLAPSE.out.gtf, ch_genome)
            ch_routea_fa = LONGREAD_GENOME_SEQ.out.fasta
        }
        ch_routea = ch_routea_fa.map { fa -> tuple('routeA', fa) }
        if (audit) {
            // residual audit: read-level removal already happened, so this is the QC gate
            AUDIT_LONGREAD(ch_routea, ch_genome, ch_contaminant, AUDIT_KEEP)
            ch_counts = ch_counts.mix(AUDIT_LONGREAD.out.counts)
            ch_routea = AUDIT_LONGREAD.out.fasta
            ch_routea_ids = AUDIT_LONGREAD.out.host_ids.first()
        }
        // every ladder rung is materialised as <label>.fasta so the reference set is a file,
        // not something the reader has to reconstruct from a work directory
        ch_base = MERGE_ISOSEQ(ch_routea.map { l, fa -> tuple('isoseq', [fa]) }).fasta
        ch_sets = ch_sets.mix(ch_base)
    }

    // ------------------------------------------ Route C: genome-guided short-read assembly
    ch_isoseq_stringtie = ch_base
    if (!params.skip_stringtie) {
        STAR_ALIGN(ch_reads, ch_star_index)
        STRINGTIE(STAR_ALIGN.out.bam)
        ch_qc = ch_qc.mix(STAR_ALIGN.out.log)

        STRINGTIE_MERGE(STRINGTIE.out.gtf.collect(), ch_collapsed_gff, ch_genome, ch_routea_ids)
        ch_st_novel = STRINGTIE_MERGE.out.fasta

        // with long reads present the short-read set only contributes the loci they missed;
        // without them it IS the base set
        ch_boosted_parts = do_routea
            ? ch_base.map { l, fa -> fa }.combine(ch_st_novel).map { base, novel -> tuple('isoseq_stringtie', [base, novel]) }
            : ch_st_novel.map { novel -> tuple('isoseq_stringtie', [novel]) }
        ch_isoseq_stringtie = MERGE_ISOSEQ_STRINGTIE(ch_boosted_parts).fasta
        ch_sets = ch_sets.mix(ch_isoseq_stringtie)
    }

    // Structural QC of the long-read models against an independent short-read annotation.
    // Reports only — the categories say which collapse products a second assembly supports.
    ch_sqanti = Channel.empty()
    if (do_routea && !params.skip_stringtie && !params.skip_sqanti) {
        SQANTI3_QC(ch_collapsed_gtf.map { gtf -> tuple('isoseq', gtf) },
                   STRINGTIE_MERGE.out.merged_gtf, ch_genome)
        ch_sqanti = SQANTI3_QC.out.summary
    }

    // -------------------------------------------------- Route D: de novo short-read assembly
    ch_denovo_novel = Channel.empty()
    if (!params.skip_denovo) {
        TRINITY(ch_reads.map { meta, reads -> reads[0] }.collect(),
                ch_reads.map { meta, reads -> reads[1] }.collect())

        ch_trin = TRINITY.out.fasta.map { fa -> tuple('denovo', fa) }
        if (audit) {
            AUDIT_DENOVO(ch_trin, ch_genome, ch_contaminant, AUDIT_KEEP)
            ch_trin = AUDIT_DENOVO.out.fasta
            ch_counts = ch_counts.mix(AUDIT_DENOVO.out.counts)
        }
        // collapse the assembler's own isoform redundancy, then keep only what the base lacks
        CLUSTER_DENOVO(ch_trin, params.cdhit_dedup_id)
        NOVEL_DENOVO(CLUSTER_DENOVO.out.fasta, ch_isoseq_stringtie, params.cdhit_dedup_id, params.prefix_denovo)
        ch_denovo_novel = NOVEL_DENOVO.out.fasta.map { l, fa -> fa }
    }

    // ------------------------------------------- Route B: reference-free long-read collapse
    ch_lrfree_novel_vs_base = Channel.empty()
    ch_lrfree_unique        = Channel.empty()
    if (do_routeb) {
        CLUSTER_LRFREE(ch_longreads.map { l, fa -> tuple('lrfree', fa) }, params.cdhit_longread_id)
        ch_lrfree = CLUSTER_LRFREE.out.fasta
        if (audit) {
            AUDIT_LRFREE(ch_lrfree, ch_genome, ch_contaminant, AUDIT_KEEP)
            ch_lrfree = AUDIT_LRFREE.out.fasta
            ch_counts = ch_counts.mix(AUDIT_LRFREE.out.counts)
        }
        CLUSTER_LRFREE_NR(ch_lrfree, params.cdhit_dedup_id)
        ch_lrfree_nr = CLUSTER_LRFREE_NR.out.fasta

        // (a) what Route B adds to the base, on its own
        NOVEL_LRFREE_A(ch_lrfree_nr, ch_isoseq_stringtie, params.cdhit_dedup_id, params.prefix_lrfree)
        ch_lrfree_novel_vs_base = NOVEL_LRFREE_A.out.fasta.map { l, fa -> fa }
    }

    // ---------------------------------------------------------------- ladder assembly
    ch_with_trinity = Channel.empty()
    if (!params.skip_denovo) {
        ch_with_trinity = MERGE_WITH_TRINITY(
            ch_isoseq_stringtie.map { l, fa -> fa }.combine(ch_denovo_novel)
                      .map { base, novel -> tuple('isoseq_stringtie_trinity', [base, novel]) }).fasta
        ch_sets = ch_sets.mix(ch_with_trinity)
    }

    ch_with_lrfree = Channel.empty()
    if (do_routeb) {
        ch_with_lrfree = MERGE_WITH_LRFREE(
            ch_isoseq_stringtie.map { l, fa -> fa }.combine(ch_lrfree_novel_vs_base)
                      .map { base, novel -> tuple('isoseq_stringtie_lrfree', [base, novel]) }).fasta
        ch_sets = ch_sets.mix(ch_with_lrfree)
    }

    // (b) what Route B still adds once the de novo set is already in — its unique contribution
    ch_all_sources = Channel.empty()
    if (do_routeb && !params.skip_denovo) {
        NOVEL_LRFREE_B(CLUSTER_LRFREE_NR.out.fasta, ch_with_trinity, params.cdhit_dedup_id, params.prefix_lrfree)
        ch_all_sources = MERGE_ALL_SOURCES(
            ch_with_trinity.map { l, fa -> fa }
                    .combine(NOVEL_LRFREE_B.out.fasta.map { l, fa -> fa })
                    .map { base, novel -> tuple('isoseq_stringtie_trinity_lrfree', [base, novel]) }).fasta
        ch_sets = ch_sets.mix(ch_all_sources)
    }

    // the widest set available becomes the input to the cross-source deduplication
    def widest_label = params.skip_denovo ? (do_routeb ? 'isoseq_stringtie_lrfree' : 'isoseq_stringtie')
                                          : (do_routeb ? 'isoseq_stringtie_trinity_lrfree'       : 'isoseq_stringtie_trinity')
    ch_widest = ch_sets.filter { l, fa -> l == widest_label }

    CLUSTER_LEAN(ch_widest.map { l, fa -> tuple('dedup', fa) }, params.cdhit_dedup_id)
    ch_lean = MERGE_DEDUP(CLUSTER_LEAN.out.fasta.map { l, fa -> tuple('all_sources_dedup', [fa]) }).fasta

    // ---------------------------------------------------------------- maps, stats, BUSCO
    // The lean set's map is the widest set's map restricted to the surviving transcripts, so
    // gene identity stays consistent between the two.
    MAKE_TX2GENE(ch_sets)
    SUBSET_TX2GENE(ch_lean, MAKE_TX2GENE.out.tx2gene.filter { l, t -> l == widest_label }.map { l, t -> t })

    ch_all_sets = ch_sets.mix(ch_lean)
    ch_tx2gene  = MAKE_TX2GENE.out.tx2gene.mix(SUBSET_TX2GENE.out.tx2gene)

    SEQ_STATS(ch_all_sets)
    ch_busco     = Channel.empty()
    ch_busco_dir = Channel.empty()
    if (!params.skip_busco) {
        BUSCO(ch_all_sets)
        ch_busco     = BUSCO.out.summary
        ch_busco_dir = BUSCO.out.dir
    }

    ROUTE_CONTRIBUTION(
        SEQ_STATS.out.stats.map { l, f -> f }.collect(),
        ch_tx2gene.map { l, f -> f }.collect(),
        params.skip_busco ? Channel.fromPath("${projectDir}/assets/NO_FILE_BUSCO").collect()
                          : ch_busco.map { l, f -> f }.collect())

    // ------------------------------------------------------- pick the quantification reference
    ch_selected = ch_all_sets.filter { l, fa -> l == params.dge_reference }
    ch_selected_t2g = ch_tx2gene.filter { l, t -> l == params.dge_reference }

    emit:
    sets       = ch_all_sets
    tx2gene    = ch_tx2gene
    selected   = ch_selected
    selected_tx2gene = ch_selected_t2g
    collapsed_gtf = ch_collapsed_gtf
    ladder     = ROUTE_CONTRIBUTION.out.table
    busco      = ch_busco
    busco_dir  = ch_busco_dir
    sqanti     = ch_sqanti
    counts     = ch_counts
    qc         = ch_qc
}
