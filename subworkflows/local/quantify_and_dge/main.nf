/*
 * Stage 5+6 — quantification, differential expression, enrichment, and the optional
 * coordinate-based rebuild of gene identity.
 *
 * The locus rebuild re-uses the SAME transcript FASTA and the SAME quantifications: only the
 * transcript->gene map changes, so re-running DGE on it is a cheap re-aggregation rather than a
 * second analysis. Both versions are kept side by side because neither map is unambiguously
 * right — see locus_rebuild_report.txt.
 */

include { SALMON_INDEX                                 } from '../../../modules/nf-core/salmon/index'
include { SALMON_INDEX as SALMON_INDEX_EXPR            } from '../../../modules/nf-core/salmon/index'
include { SALMON_QUANT                                 } from '../../../modules/nf-core/salmon/quant'
include { SALMON_QUANT as SALMON_QUANT_EXPR            } from '../../../modules/nf-core/salmon/quant'
include { SALMON_LIBTYPE_CHECK                         } from '../../../modules/local/salmon_libtype_check'
include { SALMON_LIBTYPE_CHECK as SALMON_LIBTYPE_CHECK_EXPR } from '../../../modules/local/salmon_libtype_check'
include { QUANT_SUMMARY                                } from '../../../modules/local/quant_summary'
include { FILTER_BY_EXPRESSION                         } from '../../../modules/local/filter_by_expression'
include { CORSET                                       } from '../../../modules/local/corset'
include { REMAP_ANNOTATION                             } from '../../../modules/local/remap_annotation'
include { SEQKIT_STATS as SEQ_STATS_EXPR               } from '../../../modules/nf-core/seqkit/stats'
include { TRANSCRIPTOME_BUSCO as BUSCO_EXPR            } from '../transcriptome_busco'
include { DESEQ2                                       } from '../../../modules/local/deseq2'
include { DESEQ2 as DESEQ2_LOCUS                       } from '../../../modules/local/deseq2'
include { DESEQ2 as DESEQ2_CORSET                      } from '../../../modules/local/deseq2'
include { ENRICHMENT                                   } from '../../../modules/local/enrichment'
include { ENRICHMENT as ENRICHMENT_LOCUS               } from '../../../modules/local/enrichment'
include { ENRICHMENT as ENRICHMENT_CORSET              } from '../../../modules/local/enrichment'
include { ISOFORM_SWITCHING                            } from '../../../modules/local/isoform_switching'
include { LOCUS_TX2GENE                                } from '../../../modules/local/locus_tx2gene'
include { runLocusRebuild                              } from '../utils_nfcore_hybridmetatranscriptomics_pipeline'

workflow QUANTIFY_AND_DGE {

    take:
    ch_reads          // [ meta, [r1, r2] ]
    ch_fasta          // [ meta, fasta ]
    ch_tx2gene        // [ meta, tx2gene ]
    ch_annotation     // [ meta, annotation tsv ] or empty
    ch_genome         // path (value)
    ch_collapsed_gtf  // [ meta, gtf ] or empty
    ch_samplesheet    // path (value)
    ch_contrasts      // path (value)

    main:
    ch_qc = channel.empty()
    // no genome = no decoys: the index is the transcripts alone
    // wrapped in a list: combine() spreads a list item, so a bare [] would vanish from the tuple
    ch_decoys = ch_genome.map { g -> [ g.name == 'NO_FILE_GENOME' ? [] : g ] }
    def no_ann = [ [id: 'none'], file("${projectDir}/assets/NO_FILE_ANN", checkIfExists: true) ]

    SALMON_INDEX(ch_fasta.combine(ch_decoys).map { meta, fa, g -> [ meta, fa, g ] })
    SALMON_QUANT(ch_reads, SALMON_INDEX.out.index.join(ch_tx2gene).map { meta, idx, t2g -> [ meta, idx, t2g, [] ] }.first())
    SALMON_LIBTYPE_CHECK(SALMON_QUANT.out.results, params.strandedness, params.strandedness_check)
    ch_quant_dirs = SALMON_LIBTYPE_CHECK.out.quant
    ch_quants = ch_quant_dirs.map { _meta, dir -> dir }.collect()

    /*
     * Expression filter, then a second pass. The reference the ladder produces is a union of four
     * sources with no expression evidence applied at any rung, so the first pass exists to find
     * out which transcripts the experiment actually saw; the reference, its map and its annotation
     * are cut down to those, and everything downstream is quantified against the smaller set.
     *
     * What is cut here is the DE reference, NOT the assembly. ANNOTATE and the ladder's BUSCO both
     * run upstream, on the full dge_reference set, and publish to annotation/ and transcriptome/qc/
     * — so transcripts the filter drops keep their functional annotation and still count toward the
     * reported completeness of the assembly. They are simply not quantified and not tested.
     */
    ch_ref_fasta   = ch_fasta
    ch_ref_tx2gene = ch_tx2gene
    ch_ref_ann     = ch_annotation
    ch_expr_log    = channel.empty()
    ch_expr_qc     = channel.empty()
    // no filter -> every gene in the reference is DE-eligible, signalled by the sentinel
    ch_expr_tx     = channel.value(file("${projectDir}/assets/NO_FILE_EXPRTX", checkIfExists: true))
    if (params.expr_filter_tpm > 0) {
        FILTER_BY_EXPRESSION(ch_quants, ch_fasta, ch_tx2gene, ch_annotation.ifEmpty(no_ann))
        ch_ref_fasta   = FILTER_BY_EXPRESSION.out.fasta
        ch_ref_tx2gene = FILTER_BY_EXPRESSION.out.tx2gene
        ch_expr_log    = FILTER_BY_EXPRESSION.out.log
        // DE and ORA run on the expressed genes only. With expr_filter_rescue = false (the
        // default) this is the same set as the reference and the restriction is a no-op; it stays
        // wired so that turning the rescue back on cannot silently start testing untested genes.
        ch_expr_tx     = FILTER_BY_EXPRESSION.out.expressed_tx.first()
        if (!params.skip_annotation) ch_ref_ann = FILTER_BY_EXPRESSION.out.annotation

        SALMON_INDEX_EXPR(ch_ref_fasta.combine(ch_decoys).map { meta, fa, g -> [ meta, fa, g ] })
        SALMON_QUANT_EXPR(ch_reads, SALMON_INDEX_EXPR.out.index.join(ch_ref_tx2gene).map { meta, idx, t2g -> [ meta, idx, t2g, [] ] }.first())
        SALMON_LIBTYPE_CHECK_EXPR(SALMON_QUANT_EXPR.out.results, params.strandedness, params.strandedness_check)
        ch_quant_dirs = SALMON_LIBTYPE_CHECK_EXPR.out.quant
        ch_quants = ch_quant_dirs.map { _meta, dir -> dir }.collect()

        // the gate for this rung: duplication should fall without completeness following it down
        SEQ_STATS_EXPR(ch_ref_fasta)
        ch_expr_qc = SEQ_STATS_EXPR.out.stats.map { _meta, f -> f }
        if (!params.skip_busco) {
            BUSCO_EXPR(ch_ref_fasta)
            ch_expr_qc = ch_expr_qc.mix(BUSCO_EXPR.out.summary.map { _meta, f -> f })
        }
    }

    // Salmon's own logs, so MultiQC has a quantification section
    ch_qc = ch_qc.mix(ch_quant_dirs.map { _meta, dir -> dir })

    QUANT_SUMMARY(ch_quants)

    // the gene maps themselves, so the dashboard can compare what each one calls a gene
    ch_gene_map_files = ch_ref_tx2gene.map { _meta, f -> f }

    // a third gene map from read sharing, for comparison against the other two
    ch_corset     = channel.empty()
    ch_dge_corset = channel.empty()
    ch_enr_corset = channel.empty()
    if (!params.skip_corset) {
        CORSET(ch_quants, ch_ref_tx2gene, ch_samplesheet)
        ch_corset = CORSET.out.summary
        ch_gene_map_files = ch_gene_map_files.mix(CORSET.out.tx2gene.map { _meta, f -> f })

        // Same quantifications, same code path as the other two maps — only the tx2gene differs,
        // so a difference in the results is a difference in gene identity and nothing else.
        if (!params.skip_dge) {
            DESEQ2_CORSET(ch_quants, CORSET.out.tx2gene, ch_samplesheet, ch_contrasts, ch_expr_tx)
            ch_dge_corset = DESEQ2_CORSET.out.results
            if (!params.skip_enrichment && !params.skip_annotation) {
                REMAP_ANNOTATION(CORSET.out.tx2gene, ch_ref_ann)
                ENRICHMENT_CORSET(ch_dge_corset, REMAP_ANNOTATION.out.annotation)
                ch_enr_corset = ENRICHMENT_CORSET.out.results
            }
        }
    }

    ch_dge = channel.empty()
    ch_enr = channel.empty()
    if (!params.skip_dge) {
        DESEQ2(ch_quants, ch_ref_tx2gene, ch_samplesheet, ch_contrasts, ch_expr_tx)
        ch_dge = DESEQ2.out.results
        if (!params.skip_enrichment) {
            ENRICHMENT(ch_dge, ch_ref_ann)
            ch_enr = ENRICHMENT.out.results
        }
    }

    if (!params.skip_dtu) {
        ISOFORM_SWITCHING(ch_quants, ch_ref_fasta, ch_collapsed_gtf, ch_samplesheet)
    }

    // -------------------------------------------------- gene identity from genomic loci
    ch_locus_report = channel.empty()
    ch_dge_locus    = channel.empty()
    ch_enr_locus    = channel.empty()
    if (runLocusRebuild()) {
        LOCUS_TX2GENE(ch_ref_fasta, ch_genome, ch_ref_tx2gene, ch_ref_ann)
        ch_locus_report = LOCUS_TX2GENE.out.report
        ch_gene_map_files = ch_gene_map_files.mix(LOCUS_TX2GENE.out.tx2gene.map { _meta, f -> f })

        if (!params.skip_dge) {
            DESEQ2_LOCUS(ch_quants, LOCUS_TX2GENE.out.tx2gene, ch_samplesheet, ch_contrasts, ch_expr_tx)
            ch_dge_locus = DESEQ2_LOCUS.out.results
            if (!params.skip_enrichment) {
                ENRICHMENT_LOCUS(ch_dge_locus, LOCUS_TX2GENE.out.annotation)
                ch_enr_locus = ENRICHMENT_LOCUS.out.results
            }
        }
    }

    emit:
    quants         = ch_quants        // [ quant dirs ]
    reference      = ch_ref_fasta
    expr_filter    = ch_expr_log
    expr_qc        = ch_expr_qc
    corset         = ch_corset
    gene_map_files = ch_gene_map_files
    dge_corset     = ch_dge_corset
    enrich_corset  = ch_enr_corset
    mapping        = QUANT_SUMMARY.out.table
    dge            = ch_dge
    enrichment     = ch_enr
    dge_locus      = ch_dge_locus
    enrich_locus   = ch_enr_locus
    locus_report   = ch_locus_report
    qc             = ch_qc
}
