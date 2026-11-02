/*
 * Stage 5+6 — quantification, differential expression, enrichment, and the optional
 * coordinate-based rebuild of gene identity.
 *
 * The reference this receives is already the analysed set: when the expression filter is on
 * (expr_filter_tpm > 0) EXPRESSION_FILTER has cut the ladder tip to what the experiment saw and
 * ANNOTATE has run on that subset, so there is exactly one salmon pass here and every number it
 * produces is a reported one.
 *
 * Each gene map that leaves this subworkflow also leaves it as a FASTA (GENE_FASTA), one
 * representative sequence per gene, because nothing else in the pipeline turns a gene into a
 * sequence.
 *
 * The locus rebuild re-uses the SAME transcript FASTA and the SAME quantifications: only the
 * transcript->gene map changes, so re-running DGE on it is a cheap re-aggregation rather than a
 * second analysis. Both versions are kept side by side because neither map is unambiguously
 * right — see locus_rebuild_report.txt.
 */

include { SALMON_INDEX                                 } from '../../../modules/nf-core/salmon/index'
include { SALMON_QUANT                                 } from '../../../modules/nf-core/salmon/quant'
include { SALMON_LIBTYPE_CHECK                         } from '../../../modules/local/salmon_libtype_check'
include { QUANT_SUMMARY                                } from '../../../modules/local/quant_summary'
include { CORSET                                       } from '../../../modules/local/corset'
include { REMAP_ANNOTATION                             } from '../../../modules/local/remap_annotation'
include { GENE_FASTA                                   } from '../../../modules/local/gene_fasta'
include { GENE_FASTA as GENE_FASTA_LOCUS               } from '../../../modules/local/gene_fasta'
include { GENE_FASTA as GENE_FASTA_CORSET              } from '../../../modules/local/gene_fasta'
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
    ch_fasta          // [ meta, fasta ]      the analysed reference (expression-filtered unless off)
    ch_tx2gene        // [ meta, tx2gene ]    its gene map
    ch_annotation     // [ meta, annotation tsv ] or empty
    ch_pep            // path (value)         TransDecoder peptides, or the NO_FILE sentinel
    ch_expr_tx        // path (value)         DE-eligible transcript ids, or the NO_FILE sentinel
    ch_genome         // path (value)
    ch_collapsed_gtf  // [ meta, gtf ] or empty
    ch_samplesheet    // path (value)
    ch_contrasts      // path (value)

    main:
    ch_qc = channel.empty()
    // no genome = no decoys: the index is the transcripts alone
    // wrapped in a list: combine() spreads a list item, so a bare [] would vanish from the tuple
    ch_decoys = ch_genome.map { g -> [ g.name == 'NO_FILE_GENOME' ? [] : g ] }

    SALMON_INDEX(ch_fasta.combine(ch_decoys).map { meta, fa, g -> [ meta, fa, g ] })
    SALMON_QUANT(ch_reads, SALMON_INDEX.out.index.join(ch_tx2gene).map { meta, idx, t2g -> [ meta, idx, t2g, [] ] }.first())
    SALMON_LIBTYPE_CHECK(SALMON_QUANT.out.results, params.strandedness, params.strandedness_check)
    ch_quant_dirs = SALMON_LIBTYPE_CHECK.out.quant
    ch_quants = ch_quant_dirs.map { _meta, dir -> dir }.collect()

    // Salmon's own logs, so MultiQC has a quantification section
    ch_qc = ch_qc.mix(ch_quant_dirs.map { _meta, dir -> dir })

    QUANT_SUMMARY(ch_quants)

    // the gene maps themselves, so the dashboard can compare what each one calls a gene
    ch_gene_map_files = ch_tx2gene.map { _meta, f -> f }

    // ...and the same genes as sequence. A map alone leaves a gene as a sum of counts, with no
    // FASTA anyone can BLAST, score in protein mode or attach to a paper.
    GENE_FASTA(ch_tx2gene, ch_fasta, ch_pep)
    ch_gene_rep = GENE_FASTA.out.table.map { _meta, f -> f }

    // a third gene map from read sharing, for comparison against the other two
    ch_corset     = channel.empty()
    ch_dge_corset = channel.empty()
    ch_enr_corset = channel.empty()
    if (!params.skip_corset) {
        CORSET(ch_quants, ch_tx2gene, ch_samplesheet)
        ch_corset = CORSET.out.summary
        ch_gene_map_files = ch_gene_map_files.mix(CORSET.out.tx2gene.map { _meta, f -> f })

        GENE_FASTA_CORSET(CORSET.out.tx2gene, ch_fasta, ch_pep)
        ch_gene_rep = ch_gene_rep.mix(GENE_FASTA_CORSET.out.table.map { _meta, f -> f })

        // Same quantifications, same code path as the other two maps — only the tx2gene differs,
        // so a difference in the results is a difference in gene identity and nothing else.
        if (!params.skip_dge) {
            DESEQ2_CORSET(ch_quants, CORSET.out.tx2gene, ch_samplesheet, ch_contrasts, ch_expr_tx)
            ch_dge_corset = DESEQ2_CORSET.out.results
            if (!params.skip_enrichment && !params.skip_annotation) {
                REMAP_ANNOTATION(CORSET.out.tx2gene, ch_annotation)
                ENRICHMENT_CORSET(ch_dge_corset, REMAP_ANNOTATION.out.annotation)
                ch_enr_corset = ENRICHMENT_CORSET.out.results
            }
        }
    }

    ch_dge = channel.empty()
    ch_enr = channel.empty()
    if (!params.skip_dge) {
        DESEQ2(ch_quants, ch_tx2gene, ch_samplesheet, ch_contrasts, ch_expr_tx)
        ch_dge = DESEQ2.out.results
        if (!params.skip_enrichment) {
            ENRICHMENT(ch_dge, ch_annotation)
            ch_enr = ENRICHMENT.out.results
        }
    }

    if (!params.skip_dtu) {
        ISOFORM_SWITCHING(ch_quants, ch_fasta, ch_collapsed_gtf, ch_samplesheet)
    }

    // -------------------------------------------------- gene identity from genomic loci
    ch_locus_report = channel.empty()
    ch_dge_locus    = channel.empty()
    ch_enr_locus    = channel.empty()
    if (runLocusRebuild()) {
        LOCUS_TX2GENE(ch_fasta, ch_genome, ch_tx2gene, ch_annotation)
        ch_locus_report = LOCUS_TX2GENE.out.report
        ch_gene_map_files = ch_gene_map_files.mix(LOCUS_TX2GENE.out.tx2gene.map { _meta, f -> f })

        // the primary gene definition of this pipeline, so its gene-level FASTA is the one most
        // likely to be used outside it
        GENE_FASTA_LOCUS(LOCUS_TX2GENE.out.tx2gene, ch_fasta, ch_pep)
        ch_gene_rep = ch_gene_rep.mix(GENE_FASTA_LOCUS.out.table.map { _meta, f -> f })

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
    reference      = ch_fasta
    corset         = ch_corset
    gene_map_files = ch_gene_map_files
    gene_rep       = ch_gene_rep      // one row per gene: its representative transcript and ORF
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
