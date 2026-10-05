/*
 * Stage 4 — functional annotation of the selected reference.
 *
 * The chain is ORFs -> homology (DIAMOND, eggNOG) -> domains (InterProScan) -> one table.
 * InterProScan is by far the longest step and is optional; when it is skipped the table is
 * identical minus the InterPro column, so downstream steps do not change shape.
 */

include { TRANSDECODER_LONGORF } from '../../../modules/nf-core/transdecoder/longorf'
include { TRANSDECODER_PREDICT } from '../../../modules/nf-core/transdecoder/predict'
include { CLEAN_PEP            } from '../../../modules/local/clean_pep'
include { DIAMOND_BLASTP       } from '../../../modules/nf-core/diamond/blastp'
include { EGGNOGMAPPER         } from '../../../modules/nf-core/eggnogmapper'
include { INTERPROSCAN         } from '../../../modules/nf-core/interproscan'
include { INTEGRATE_ANNOTATION } from '../../../modules/local/integrate_annotation'
include { isUnset              } from '../utils_nfcore_hybridmetatranscriptomics_pipeline'

workflow ANNOTATE {

    take:
    ch_fasta      // [ meta, fasta ]
    ch_tx2gene    // [ meta, tx2gene ]
    ch_dbs        // [ meta, dmnd ]

    main:
    // TransDecoder resumes from checkpoints purely on file presence, so a stale working dir
    // silently emits peptides keyed on a previous reference's ids. Nextflow gives each task a
    // clean directory, which removes that failure mode by construction.
    TRANSDECODER_LONGORF(ch_fasta)
    TRANSDECODER_PREDICT(ch_fasta, TRANSDECODER_LONGORF.out.folder)
    CLEAN_PEP(TRANSDECODER_PREDICT.out.pep)
    ch_pep = CLEAN_PEP.out.pep

    ch_search = ch_pep.combine(ch_dbs).multiMap { meta, pep, meta2, db ->
        query: [ meta, pep ]
        db:    [ meta2, db ]
    }
    DIAMOND_BLASTP(ch_search.query, ch_search.db, 6, 'qseqid sseqid pident length evalue bitscore stitle')
    ch_diamond = DIAMOND_BLASTP.out.txt.map { _meta, txt -> txt }.collect()

    ch_emapper = channel.value(file("${projectDir}/assets/NO_FILE", checkIfExists: true))
    if (!params.skip_eggnog) {
        EGGNOGMAPPER(ch_pep, [ 'diamond', [] ], file(params.eggnog_data_dir, checkIfExists: true))
        ch_emapper = EGGNOGMAPPER.out.annotations.map { _meta, f -> f }.first()
    }

    ch_interpro = channel.value(file("${projectDir}/assets/NO_FILE_IPR", checkIfExists: true))
    if (!params.skip_interproscan) {
        INTERPROSCAN(ch_pep, isUnset(params.interproscan_db) ? [] : file(params.interproscan_db, checkIfExists: true))
        ch_interpro = INTERPROSCAN.out.tsv.map { _meta, f -> f }.first()
    }

    // join on the reference so the table cannot be built from a mismatched pair
    INTEGRATE_ANNOTATION(
        ch_tx2gene.join(TRANSDECODER_PREDICT.out.pep),
        ch_diamond,
        ch_emapper,
        ch_interpro)

    emit:
    table   = INTEGRATE_ANNOTATION.out.table       // [ meta, tsv ]
    summary = INTEGRATE_ANNOTATION.out.summary
    pep     = TRANSDECODER_PREDICT.out.pep
    cds     = TRANSDECODER_PREDICT.out.cds
}
