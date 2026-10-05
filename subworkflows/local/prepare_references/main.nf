/*
 * Stage 0 — everything the run needs before it can touch a read.
 *
 * Each reference is taken from a local path when one is given, and downloaded otherwise. A new
 * site therefore needs nothing staged by hand: only the accessions/URLs in the config.
 */

include { DOWNLOAD_GENOME              } from '../../../modules/local/download_genome'
include { DOWNLOAD_CONTAMINANT_GENOMES } from '../../../modules/local/download_contaminant_genomes'
include { COMBINE_CONTAMINANTS         } from '../../../modules/local/combine_contaminants'
include { DOWNLOAD_RRNA_DB             } from '../../../modules/local/download_rrna_db'
include { DOWNLOAD_PROTEIN_DB          } from '../../../modules/local/download_protein_db'
include { GUNZIP as GUNZIP_GENOME      } from '../../../modules/nf-core/gunzip'
include { DIAMOND_MAKEDB               } from '../../../modules/nf-core/diamond/makedb'
include { STAR_GENOMEGENERATE          } from '../../../modules/nf-core/star/genomegenerate'
include { isUnset; runDecontam; runStringtie } from '../utils_nfcore_hybridmetatranscriptomics_pipeline'

workflow PREPARE_REFERENCES {

    main:

    // ---------------------------------------------------------------- host genome
    if (!isUnset(params.genome)) {
        if (params.genome.toString().endsWith('.gz')) {
            // STAR, pbmm2's index and gffread all want it uncompressed
            ch_genome = GUNZIP_GENOME(channel.value([ [id: 'genome'], file(params.genome, checkIfExists: true) ]))
                .gunzip.map { _meta, fa -> fa }.first()
        } else {
            ch_genome = channel.value(file(params.genome, checkIfExists: true))
        }
    } else if (!isUnset(params.genome_accession)) {
        ch_genome = DOWNLOAD_GENOME(params.genome_accession, params.genome_name).fasta.first()
    } else {
        // No genome: the genome-guided steps are off and the steps that can run without one see
        // this empty placeholder — decontamination then removes whatever places on a
        // contaminant, and salmon indexes without decoys.
        ch_genome = channel.value(file("${projectDir}/assets/NO_FILE_GENOME", checkIfExists: true))
    }

    // ------------------------------------------------------- contaminant genomes
    // Any one of: a combined FASTA, a directory of FASTAs, or a list of NCBI accessions.
    // None of them = no decontamination; with --skip_decontam they are not fetched either.
    if (!runDecontam()) {
        ch_contam = channel.empty()
    } else if (!isUnset(params.contaminant_fasta)) {
        ch_contam = channel.value(file(params.contaminant_fasta, checkIfExists: true))
    } else if (!isUnset(params.contaminant_dir)) {
        ch_contam = COMBINE_CONTAMINANTS(file(params.contaminant_dir, checkIfExists: true)).fasta.first()
    } else {
        ch_contam = DOWNLOAD_CONTAMINANT_GENOMES(params.contaminant_accessions).fasta.first()
    }

    // ------------------------------------------------------------- rRNA database
    if (params.skip_rrna) {
        ch_rrna = channel.empty()
    } else if (!isUnset(params.rrna_fasta)) {
        ch_rrna = channel.value(file(params.rrna_fasta, checkIfExists: true))
    } else {
        ch_rrna = DOWNLOAD_RRNA_DB(params.rrna_db_base, params.rrna_db_sets).fasta.first()
    }

    // --------------------------------------------------------------- STAR index
    // No GTF: junctions are discovered de novo by the 2-pass alignment. An unreliable annotation
    // of a draft genome would bias the discovery this pipeline exists to perform.
    if (!runStringtie()) {
        ch_star = channel.empty()
    } else if (!isUnset(params.star_index)) {
        ch_star = channel.value([ [id: 'genome'], file(params.star_index, checkIfExists: true) ])
    } else {
        ch_star = STAR_GENOMEGENERATE(ch_genome.map { fa -> [ [id: 'genome'], fa ] }, [ [:], [] ]).index.first()
    }

    // ------------------------------------------------------ protein DBs (DIAMOND)
    if (params.skip_annotation) {
        ch_dbs = channel.empty()
    } else if (!isUnset(params.diamond_db_dir)) {
        ch_dbs = channel.fromPath("${params.diamond_db_dir}/*.dmnd", checkIfExists: true)
            .map { db -> [ [id: db.simpleName], db ] }
    } else {
        def urls = []
        if (!isUnset(params.swissprot_url))  urls << [ [id: 'uniprot_sprot'], params.swissprot_url ]
        if (!isUnset(params.specialist_url)) urls << [ [id: params.specialist_name], params.specialist_url ]
        DOWNLOAD_PROTEIN_DB(channel.fromList(urls))
        ch_dbs = DIAMOND_MAKEDB(DOWNLOAD_PROTEIN_DB.out.fasta, [], [], []).db
    }

    emit:
    genome      = ch_genome      // path (value): FASTA or the NO_FILE_GENOME placeholder
    contaminant = ch_contam      // path (value)
    rrna        = ch_rrna        // path (value)
    star_index  = ch_star        // [ meta, dir ] (value)
    diamond_dbs = ch_dbs         // [ meta, dmnd ]
}
