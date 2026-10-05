/*
 * Stage 1+2 (short reads) — from raw reads to host-only reads.
 *
 *   merge lanes -> fastp -> rRNA depletion -> FastQC -> contaminant removal
 *
 * Contaminant removal happens at the READ level, before anything is assembled: an assembly built
 * from mixed reads cannot be cleanly separated afterwards, because chimeric contigs and shared
 * conserved regions have already been baked in.
 */

include { CAT_FASTQ      } from '../../../modules/nf-core/cat/fastq'
include { FASTP          } from '../../../modules/nf-core/fastp'
include { BBMAP_BBDUK    } from '../../../modules/nf-core/bbmap/bbduk'
include { FASTQC         } from '../../../modules/nf-core/fastqc'
include { RRNA_SUMMARY   } from '../../../modules/local/rrna_summary'
include { DECONTAM_SHORT } from '../../../modules/local/decontam_short'
include { runDecontam    } from '../utils_nfcore_hybridmetatranscriptomics_pipeline'

workflow SHORT_READ_PREPROCESS {

    take:
    ch_units          // [ meta, [ r1_L1, r2_L1, r1_L2, r2_L2, ... ] ]
    ch_rrna           // path (value)
    ch_decon_index    // path (value), or empty when decontamination is off

    main:
    ch_qc  = channel.empty()
    ch_tsv = channel.empty()

    // one row per lane in the samplesheet; gzip members concatenate legally
    CAT_FASTQ(ch_units)
    ch_reads = CAT_FASTQ.out.reads

    if (!params.skip_fastp) {
        FASTP(ch_reads.map { meta, reads -> [ meta, reads, [] ] }, false, false, false)
        ch_reads = FASTP.out.reads
        ch_qc    = ch_qc.mix(FASTP.out.json.map { _meta, json -> json })
    }

    if (!params.skip_rrna) {
        BBMAP_BBDUK(ch_reads, ch_rrna)
        ch_reads = BBMAP_BBDUK.out.reads
        // per-sample rRNA percentage, parsed from bbduk's own accounting
        RRNA_SUMMARY(BBMAP_BBDUK.out.log)
        ch_qc  = ch_qc.mix(BBMAP_BBDUK.out.log.map { _meta, log -> log })
        ch_tsv = ch_tsv.mix(RRNA_SUMMARY.out.summary.collect().map { t -> [ 'rrna_per_sample', t ] })
    }

    FASTQC(ch_reads)
    ch_qc = ch_qc.mix(FASTQC.out.zip.map { _meta, zip -> zip })

    if (runDecontam()) {
        DECONTAM_SHORT(ch_reads, ch_decon_index, params.contaminant_prefix)
        ch_reads = DECONTAM_SHORT.out.reads
        ch_tsv   = ch_tsv.mix(DECONTAM_SHORT.out.summary.collect().map { t -> [ 'decontam_per_sample', t ] })
    }

    emit:
    reads  = ch_reads    // [ meta, [ r1, r2 ] ]
    qc     = ch_qc       // files for MultiQC
    tables = ch_tsv      // [ table_name, [ per-sample tsv ] ]
}
