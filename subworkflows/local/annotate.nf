/*
 * Stage 4 — functional annotation of the selected reference.
 *
 * The chain is ORFs -> homology (DIAMOND, eggNOG) -> domains (InterProScan) -> one table.
 * InterProScan is by far the longest step and is optional; when it is skipped the table is
 * identical minus the InterPro column, so downstream steps do not change shape.
 */

include { TRANSDECODER          } from '../../modules/local/annotation'
include { CLEAN_PEP             } from '../../modules/local/annotation'
include { DIAMOND_BLASTP        } from '../../modules/local/annotation'
include { EGGNOG_MAPPER         } from '../../modules/local/annotation'
include { INTERPROSCAN          } from '../../modules/local/annotation'
include { INTEGRATE_ANNOTATION  } from '../../modules/local/annotation'

workflow ANNOTATE {

    take:
    ch_fasta      // [ label, fasta ]
    ch_tx2gene    // [ label, tx2gene ]
    ch_dbs        // [ db_name, dmnd ]

    main:
    TRANSDECODER(ch_fasta)
    CLEAN_PEP(TRANSDECODER.out.pep)
    ch_pep = CLEAN_PEP.out.pep

    DIAMOND_BLASTP(ch_pep.combine(ch_dbs))
    ch_diamond = DIAMOND_BLASTP.out.hits.map { name, tsv -> tsv }.collect()

    ch_emapper = Channel.fromPath("${projectDir}/assets/NO_FILE")
    if (!params.skip_eggnog) {
        if (!params.eggnog_data_dir)
            error "eggNOG-mapper needs --eggnog_data_dir (or run with --skip_eggnog)."
        ch_emapper = EGGNOG_MAPPER(ch_pep, params.eggnog_data_dir).annotations
    }

    ch_interpro = Channel.fromPath("${projectDir}/assets/NO_FILE_IPR")
    if (!params.skip_interproscan) {
        ch_interpro = INTERPROSCAN(ch_pep).tsv
    }

    // join on the reference label so the table cannot be built from a mismatched pair
    INTEGRATE_ANNOTATION(
        ch_tx2gene.join(TRANSDECODER.out.pep),
        ch_diamond,
        ch_emapper.first(),
        ch_interpro.first())

    emit:
    table   = INTEGRATE_ANNOTATION.out.table       // [ label, tsv ]
    summary = INTEGRATE_ANNOTATION.out.summary
    pep     = TRANSDECODER.out.pep
    cds     = TRANSDECODER.out.cds
}
