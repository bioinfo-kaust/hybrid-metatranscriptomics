/*
 * Stage 3b — the expression filter that defines the analysed reference.
 *
 * The ladder applies no expression evidence at any rung, so its tip is a cross-source union that
 * carries every near-duplicate and every unexpressed fragment the four assembly routes produced —
 * which is what 31% BUSCO duplication looks like. One salmon pass over that tip measures what the
 * experiment actually saw; the reference and its gene map are cut to the survivors, and everything
 * after this point — annotation, quantification, DE, enrichment, figures and the reports — runs on
 * that one set.
 *
 * Nothing here is a reported number. This pass exists only to produce TPMs for the filter: its
 * quantifications are not published, and the counts every DE number comes from are produced by the
 * single salmon pass in QUANTIFY_AND_DGE, against the filtered reference.
 *
 * What the filter cuts is the ANALYSED set, not the assembly. BUSCO and seqkit still run on every
 * ladder rung in BUILD_TRANSCRIPTOME, so the completeness and size quoted for the assembly remain
 * properties of the full deduplicated set; this subworkflow adds the same two measurements for the
 * filtered set, so the two can be read on one line each.
 */

include { SALMON_INDEX        as SALMON_INDEX_PRE        } from '../../../modules/nf-core/salmon/index'
include { SALMON_QUANT        as SALMON_QUANT_PRE        } from '../../../modules/nf-core/salmon/quant'
include { SALMON_LIBTYPE_CHECK as SALMON_LIBTYPE_CHECK_PRE } from '../../../modules/local/salmon_libtype_check'
include { FILTER_BY_EXPRESSION                          } from '../../../modules/local/filter_by_expression'
include { SEQKIT_STATS        as SEQ_STATS_EXPR         } from '../../../modules/nf-core/seqkit/stats'
include { TRANSCRIPTOME_BUSCO as BUSCO_EXPR             } from '../transcriptome_busco'

workflow EXPRESSION_FILTER {

    take:
    ch_reads          // [ meta, [r1, r2] ]
    ch_fasta          // [ meta, fasta ]      the selected ladder rung
    ch_tx2gene        // [ meta, tx2gene ]    its gene map
    ch_genome         // path (value)

    main:
    // no genome = no decoys: the index is the transcripts alone
    // wrapped in a list: combine() spreads a list item, so a bare [] would vanish from the tuple
    ch_decoys = ch_genome.map { g -> [ g.name == 'NO_FILE_GENOME' ? [] : g ] }

    SALMON_INDEX_PRE(ch_fasta.combine(ch_decoys).map { meta, fa, g -> [ meta, fa, g ] })
    SALMON_QUANT_PRE(ch_reads, SALMON_INDEX_PRE.out.index.join(ch_tx2gene).map { meta, idx, t2g -> [ meta, idx, t2g, [] ] }.first())
    // the same library-type check as the reported pass: a mis-detected library here would filter
    // on the wrong TPMs and the error would be invisible by the time anything is tested
    SALMON_LIBTYPE_CHECK_PRE(SALMON_QUANT_PRE.out.results, params.strandedness, params.strandedness_check)
    ch_quants = SALMON_LIBTYPE_CHECK_PRE.out.quant.map { _meta, dir -> dir }.collect()

    FILTER_BY_EXPRESSION(ch_quants, ch_fasta, ch_tx2gene)

    // the gate for this rung: duplication should fall without completeness following it down
    SEQ_STATS_EXPR(FILTER_BY_EXPRESSION.out.fasta)
    ch_qc = SEQ_STATS_EXPR.out.stats.map { _meta, f -> f }
    if (!params.skip_busco) {
        BUSCO_EXPR(FILTER_BY_EXPRESSION.out.fasta)
        ch_qc = ch_qc.mix(BUSCO_EXPR.out.summary.map { _meta, f -> f })
    }

    emit:
    fasta        = FILTER_BY_EXPRESSION.out.fasta         // [ meta, <rung>_expressed.fasta ]
    tx2gene      = FILTER_BY_EXPRESSION.out.tx2gene       // [ meta, tx2gene.<rung>_expressed.tsv ]
    expressed_tx = FILTER_BY_EXPRESSION.out.expressed_tx  // the transcript ids that passed
    filter_log   = FILTER_BY_EXPRESSION.out.log           // the filter's own accounting
    qc           = ch_qc                                  // stats + BUSCO of the filtered set
}
