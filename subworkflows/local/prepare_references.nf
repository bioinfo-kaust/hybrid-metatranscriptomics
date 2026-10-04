/*
 * Stage 0 — everything the run needs before it can touch a read.
 *
 * Each reference is taken from a local path when one is given, and downloaded otherwise. A new
 * site therefore needs nothing staged by hand: only the accessions/URLs in the config.
 */

include { DOWNLOAD_GENOME              } from '../../modules/local/references'
include { DOWNLOAD_CONTAMINANT_GENOMES } from '../../modules/local/references'
include { COMBINE_CONTAMINANTS         } from '../../modules/local/references'
include { DOWNLOAD_RRNA_DB             } from '../../modules/local/references'
include { DOWNLOAD_PROTEIN_DB          } from '../../modules/local/references'
include { MAKE_DIAMOND_DB              } from '../../modules/local/references'
include { STAR_INDEX                   } from '../../modules/local/references'

// Nextflow passes command-line values as strings, so `--contaminant_fasta null` arrives as the
// four-character string "null". Treat that (and an empty string) as "not set".
def unset(v) { (v == null || v == 'null' || v.toString().trim() == '') ? null : v }

workflow PREPARE_REFERENCES {

    main:

    // ---------------------------------------------------------------- host genome
    if (unset(params.genome)) {
        ch_genome = Channel.fromPath(params.genome, checkIfExists: true).first()
    } else if (unset(params.genome_accession)) {
        ch_genome = DOWNLOAD_GENOME(params.genome_accession, params.genome_name).fasta.first()
    } else {
        error "Provide --genome (a FASTA) or --genome_accession (downloaded with `datasets`)."
    }

    // ------------------------------------------------------- contaminant genomes
    // Any one of: a combined FASTA, a directory of FASTAs, or a list of NCBI accessions.
    if (unset(params.contaminant_fasta)) {
        ch_contam = Channel.fromPath(params.contaminant_fasta, checkIfExists: true).first()
    } else if (unset(params.contaminant_dir)) {
        ch_contam = COMBINE_CONTAMINANTS(Channel.fromPath(params.contaminant_dir, checkIfExists: true)).fasta.first()
    } else if (unset(params.contaminant_accessions)) {
        ch_contam = DOWNLOAD_CONTAMINANT_GENOMES(unset(params.contaminant_accessions)).fasta.first()
    } else {
        ch_contam = Channel.empty()
    }

    // ------------------------------------------------------------- rRNA database
    if (params.skip_rrna) {
        ch_rrna = Channel.empty()
    } else if (unset(params.rrna_fasta)) {
        ch_rrna = Channel.fromPath(params.rrna_fasta, checkIfExists: true).first()
    } else {
        ch_rrna = DOWNLOAD_RRNA_DB().fasta.first()
    }

    // --------------------------------------------------------------- STAR index
    if (params.skip_stringtie) {
        ch_star = Channel.empty()
    } else if (unset(params.star_index)) {
        ch_star = Channel.fromPath(params.star_index, checkIfExists: true).first()
    } else {
        ch_star = STAR_INDEX(ch_genome).index.first()
    }

    // ------------------------------------------------------ protein DBs (DIAMOND)
    if (params.skip_annotation) {
        ch_dbs = Channel.empty()
    } else if (unset(params.diamond_db_dir)) {
        ch_dbs = Channel.fromPath("${unset(params.diamond_db_dir)}/*.dmnd", checkIfExists: true)
                        .map { db -> tuple(db.simpleName, db) }
    } else {
        def urls = []
        if (params.swissprot_url)  urls << tuple('uniprot_sprot', params.swissprot_url)
        if (params.specialist_url) urls << tuple(params.specialist_name, params.specialist_url)
        ch_dbs = MAKE_DIAMOND_DB(DOWNLOAD_PROTEIN_DB(Channel.fromList(urls)).fasta).db
    }

    emit:
    genome      = ch_genome
    contaminant = ch_contam
    rrna        = ch_rrna
    star_index  = ch_star
    diamond_dbs = ch_dbs
}
