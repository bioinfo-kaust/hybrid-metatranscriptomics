/*
 * BUSCO completeness of one or more transcript sets.
 *
 * Headers are replaced by serial ids first: only the completeness string is consumed downstream,
 * and BUSCO crashes (KeyError in GenomeAnalysis.get_exon_records) when metaeuk truncates very long
 * de novo headers into colliding gene ids.
 */

include { SEQKIT_REPLACE as BUSCO_RENAME } from '../../../modules/nf-core/seqkit/replace'
include { BUSCO_BUSCO                    } from '../../../modules/nf-core/busco/busco'
include { BUSCO_SUMMARY                  } from '../../../modules/local/busco_summary'
include { isUnset                        } from '../utils_nfcore_hybridmetatranscriptomics_pipeline'

workflow TRANSCRIPTOME_BUSCO {

    take:
    ch_fasta          // [ meta, fasta ]

    main:
    BUSCO_RENAME(ch_fasta, 'fasta')
    BUSCO_BUSCO(
        BUSCO_RENAME.out.fastx,
        params.busco_mode,
        params.busco_lineage,
        isUnset(params.busco_download_path) ? [] : file(params.busco_download_path, checkIfExists: true),
        [],
        true)
    // every set gets a line: a BUSCO that failed (errorStrategy ignore) or wrote no summary
    // reports NA instead of silently holding up everything that waits on the summary.
    // The full ortholog table rides along on the same join, for the same reason — a rung without
    // one must not stall the rungs that have one.
    BUSCO_SUMMARY(
        ch_fasta
            .map { meta, _fa -> [ meta.id, meta ] }
            .join(BUSCO_BUSCO.out.short_summaries_txt.map { meta, txt -> [ meta.id, txt ] }, remainder: true)
            .join(BUSCO_BUSCO.out.full_table.map { meta, ft -> [ meta.id, ft ] }, remainder: true)
            .filter { _id, meta, _txt, _ft -> meta != null }
            .map { _id, meta, txt, ft -> [ meta, txt ?: [], ft ?: [] ] }
    )

    emit:
    summary = BUSCO_SUMMARY.out.summary                 // [ meta, <id>.busco.txt ]
    // which ortholog each sequence carries: the label-free truth the gene-map evaluation uses.
    // Ids are BUSCO_RENAME's serial ids, resolved through the input FASTA's own order.
    // Taken from BUSCO_SUMMARY, not BUSCO_BUSCO: see that module's header for why.
    full_table = BUSCO_SUMMARY.out.full_table           // [ meta, <id>.full_table.tsv ]
    // only the specific-lineage summaries: MultiQC parses these, and the generic ones repeat them
    mqc     = BUSCO_BUSCO.out.short_summaries_txt
        .flatMap { _meta, files -> [files].flatten() }
        .filter { f -> f.name ==~ /short_summary\.specific\..*\.txt/ }
}
