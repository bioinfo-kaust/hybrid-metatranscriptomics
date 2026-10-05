/*
 * Stage 5+6 — quantification, differential expression, enrichment, and the optional
 * coordinate-based rebuild of gene identity.
 *
 * The locus rebuild re-uses the SAME transcript FASTA and the SAME quantifications: only the
 * transcript->gene map changes, so re-running DGE on it is a cheap re-aggregation rather than a
 * second analysis. Both versions are kept side by side because neither map is unambiguously
 * right — see locus_rebuild_report.txt.
 */

include { SALMON_INDEX      } from '../../modules/local/quant'
include { SALMON_INDEX as SALMON_INDEX_EXPR } from '../../modules/local/quant'
include { SALMON_QUANT      } from '../../modules/local/quant'
include { SALMON_QUANT as SALMON_QUANT_EXPR } from '../../modules/local/quant'
include { QUANT_SUMMARY     } from '../../modules/local/quant'
include { FILTER_BY_EXPRESSION } from '../../modules/local/quant'
include { CORSET            } from '../../modules/local/quant'
include { REMAP_ANNOTATION  } from '../../modules/local/quant'
include { SEQ_STATS as SEQ_STATS_EXPR } from '../../modules/local/references'
include { BUSCO as BUSCO_EXPR         } from '../../modules/local/assembly'
include { DESEQ2            } from '../../modules/local/dge'
include { DESEQ2 as DESEQ2_LOCUS } from '../../modules/local/dge'
include { DESEQ2 as DESEQ2_CORSET } from '../../modules/local/dge'
include { ENRICHMENT        } from '../../modules/local/dge'
include { ENRICHMENT as ENRICHMENT_LOCUS } from '../../modules/local/dge'
include { ENRICHMENT as ENRICHMENT_CORSET } from '../../modules/local/dge'
include { ISOFORM_SWITCHING } from '../../modules/local/dge'
include { LOCUS_TX2GENE     } from '../../modules/local/dge'

workflow QUANTIFY_AND_DGE {

    take:
    ch_reads          // [ meta, [r1, r2] ]
    ch_fasta          // [ label, fasta ]
    ch_tx2gene        // [ label, tx2gene ]
    ch_annotation     // [ label, annotation tsv ]
    ch_genome
    ch_collapsed_gtf
    ch_samplesheet
    ch_contrasts

    main:
    SALMON_INDEX(ch_fasta, ch_genome)
    SALMON_QUANT(ch_reads, SALMON_INDEX.out.index.first())
    ch_quants = SALMON_QUANT.out.quant.map { meta, dir -> dir }.collect()

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
    ch_expr_log    = Channel.empty()
    ch_expr_qc     = Channel.empty()
    // no filter -> every gene in the reference is DE-eligible, signalled by the sentinel
    ch_expr_tx     = Channel.value(file("${projectDir}/assets/NO_FILE_EXPRTX"))
    if (params.expr_filter_tpm > 0) {
        FILTER_BY_EXPRESSION(ch_quants, ch_fasta, ch_tx2gene,
                             ch_annotation.ifEmpty { tuple('none', file("${projectDir}/assets/NO_FILE_ANN")) })
        ch_ref_fasta   = FILTER_BY_EXPRESSION.out.fasta
        ch_ref_tx2gene = FILTER_BY_EXPRESSION.out.tx2gene
        ch_expr_log    = FILTER_BY_EXPRESSION.out.log
        // DE and ORA run on the expressed genes only. With expr_filter_rescue = false (the
        // default) this is the same set as the reference and the restriction is a no-op; it stays
        // wired so that turning the rescue back on cannot silently start testing untested genes.
        ch_expr_tx     = FILTER_BY_EXPRESSION.out.expressed_tx.first()
        if (!params.skip_annotation) ch_ref_ann = FILTER_BY_EXPRESSION.out.annotation

        SALMON_INDEX_EXPR(ch_ref_fasta, ch_genome)
        SALMON_QUANT_EXPR(ch_reads, SALMON_INDEX_EXPR.out.index.first())
        ch_quants = SALMON_QUANT_EXPR.out.quant.map { meta, dir -> dir }.collect()

        // the gate for this rung: duplication should fall without completeness following it down
        SEQ_STATS_EXPR(ch_ref_fasta)
        ch_expr_qc = SEQ_STATS_EXPR.out.stats.map { l, f -> f }
        if (!params.skip_busco) {
            ch_expr_qc = ch_expr_qc.mix(BUSCO_EXPR(ch_ref_fasta).summary.map { l, f -> f })
        }
    }

    QUANT_SUMMARY(ch_quants)

    // the gene maps themselves, so the dashboard can compare what each one calls a gene
    ch_gene_map_files = ch_ref_tx2gene.map { l, f -> f }

    // a third gene map from read sharing, for comparison against the other two
    ch_corset       = Channel.empty()
    ch_dge_corset   = Channel.empty()
    ch_enr_corset   = Channel.empty()
    if (!params.skip_corset) {
        CORSET(ch_quants, ch_ref_tx2gene, ch_samplesheet)
        ch_corset = CORSET.out.summary
        ch_gene_map_files = ch_gene_map_files.mix(CORSET.out.tx2gene.map { l, f -> f })

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

    ch_dge = Channel.empty()
    ch_enr = Channel.empty()
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
    ch_locus_report = Channel.empty()
    ch_dge_locus    = Channel.empty()
    ch_enr_locus    = Channel.empty()
    if (Inputs.locusRebuild(params)) {
        LOCUS_TX2GENE(ch_ref_fasta, ch_genome, ch_ref_tx2gene, ch_ref_ann)
        ch_locus_report = LOCUS_TX2GENE.out.report
        ch_gene_map_files = ch_gene_map_files.mix(LOCUS_TX2GENE.out.tx2gene.map { l, f -> f })

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
    quants       = ch_quants
    reference    = ch_ref_fasta
    expr_filter  = ch_expr_log
    expr_qc      = ch_expr_qc
    corset       = ch_corset
    gene_map_files = ch_gene_map_files
    dge_corset   = ch_dge_corset
    enrich_corset = ch_enr_corset
    mapping      = QUANT_SUMMARY.out.table
    dge          = ch_dge
    enrichment   = ch_enr
    dge_locus    = ch_dge_locus
    enrich_locus = ch_enr_locus
    locus_report = ch_locus_report
}
