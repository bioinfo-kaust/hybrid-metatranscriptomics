/*
 * Stage 1+2 — from raw reads to host-only reads.
 *
 *   short reads:  merge units -> fastp -> rRNA depletion -> FastQC -> contaminant removal
 *   long reads :  normalise  -> rRNA depletion -> contaminant classification
 *
 * Contaminant removal happens at the READ level, before anything is assembled: an assembly built
 * from mixed reads cannot be cleanly separated afterwards, because chimeric contigs and shared
 * conserved regions have already been baked in.
 */

include { CAT_FASTQ      } from '../../modules/local/shortread_qc'
include { FASTP          } from '../../modules/local/shortread_qc'
include { BBDUK_RRNA_PE  } from '../../modules/local/shortread_qc'
include { BBDUK_RRNA_SE  } from '../../modules/local/shortread_qc'
include { FASTQC         } from '../../modules/local/shortread_qc'
include { LONGREAD_PREP  } from '../../modules/local/shortread_qc'
include { LONGREAD_MEMBERSHIP } from '../../modules/local/shortread_qc'
include { LONGREAD_POOL  } from '../../modules/local/shortread_qc'
include { BUILD_DECON_INDEX } from '../../modules/local/decontam'
include { DECONTAM_SHORT    } from '../../modules/local/decontam'
include { DECONTAM_LONG     } from '../../modules/local/decontam'

// `--longread_membership null` on the command line arrives as the four-character string "null"
def unset_lr(v) { (v == null || v == 'null' || v.toString().trim() == '') ? null : v }

workflow SHORT_READ_PREPROCESS {

    take:
    ch_units          // [ meta, [r1...], [r2...] ]
    ch_rrna
    ch_decon_index

    main:
    ch_qc  = Channel.empty()
    ch_tsv = Channel.empty()

    // one row per unit in the sample sheet = one lane; gzip members concatenate legally
    CAT_FASTQ(ch_units)
    ch_reads = CAT_FASTQ.out.reads

    if (!params.skip_fastp) {
        FASTP(ch_reads)
        ch_reads = FASTP.out.reads
        ch_qc = ch_qc.mix(FASTP.out.json)
    }

    if (!params.skip_rrna) {
        BBDUK_RRNA_PE(ch_reads, ch_rrna)
        ch_reads = BBDUK_RRNA_PE.out.reads
        ch_qc  = ch_qc.mix(BBDUK_RRNA_PE.out.stats, BBDUK_RRNA_PE.out.log)
        ch_tsv = ch_tsv.mix(BBDUK_RRNA_PE.out.summary.collect().map { t -> tuple('rrna_per_sample', t) })
    }

    FASTQC(ch_reads)
    ch_qc = ch_qc.mix(FASTQC.out.report)

    ch_decontam_tsv = Channel.empty()
    if (!params.skip_decontam) {
        DECONTAM_SHORT(ch_reads, ch_decon_index)
        ch_reads = DECONTAM_SHORT.out.reads
        ch_decontam_tsv = DECONTAM_SHORT.out.summary.collect().map { t -> tuple('decontam_per_sample', t) }
        ch_tsv = ch_tsv.mix(ch_decontam_tsv)
    }

    emit:
    reads    = ch_reads
    qc       = ch_qc
    tables   = ch_tsv
    decontam = ch_decontam_tsv
}

workflow LONG_READ_PREPROCESS {

    take:
    ch_long           // [ meta, reads ]
    ch_rrna
    ch_genome
    ch_contaminant

    main:
    ch_qc     = Channel.empty()
    ch_tables = Channel.empty()

    LONGREAD_PREP(ch_long)
    ch_seqs = LONGREAD_PREP.out.seqs

    // which sequence came from which sample, taken before anything is filtered so the per-sample
    // input counts are the real ones. An external map (--longread_membership) covers the case
    // where the long reads were already pooled before the pipeline saw them.
    def placeholder = file("${projectDir}/assets/NO_FILE_MEMBERSHIP")
    ch_membership = unset_lr(params.longread_membership)
        ? Channel.fromPath(unset_lr(params.longread_membership), checkIfExists: true).collect()
        : LONGREAD_MEMBERSHIP(ch_seqs).tsv.collect().ifEmpty([placeholder])

    ch_rrna_hits = Channel.fromPath("${projectDir}/assets/NO_FILE_RRNA").collect()
    if (!params.skip_rrna) {
        BBDUK_RRNA_SE(ch_seqs, ch_rrna)
        ch_seqs = BBDUK_RRNA_SE.out.seqs
        ch_qc = ch_qc.mix(BBDUK_RRNA_SE.out.stats)
        ch_tables = ch_tables.mix(BBDUK_RRNA_SE.out.summary.collect().map { t -> tuple('longread_rrna', t) })
        ch_rrna_hits = BBDUK_RRNA_SE.out.removed.collect()
    }

    // sample -> condition, for the per-sample table
    ch_conditions = ch_long
        .map { meta, reads -> "${meta.id}\t${meta.condition}\n" }
        .collectFile(name: 'longread_conditions.tsv', sort: true)

    // The reference is built from all long reads together: these libraries are usually few and
    // unreplicated, and isoform discovery benefits from the pooled depth.
    ch_pooled = LONGREAD_POOL(ch_seqs.map { meta, seqs -> seqs }.collect()).seqs
                              .map { fa -> tuple('longreads', fa) }

    if (!params.skip_decontam) {
        DECONTAM_LONG(ch_pooled, ch_genome, ch_contaminant, ch_membership, ch_rrna_hits, ch_conditions)
        ch_pooled = DECONTAM_LONG.out.seqs
        ch_tables = ch_tables
            .mix(DECONTAM_LONG.out.summary.map { t -> tuple('longread_decontam', [t]) })
            .mix(DECONTAM_LONG.out.per_sample.map { t -> tuple('longread_per_sample', [t]) })
    }

    emit:
    seqs   = ch_pooled          // [ 'longreads', fasta ]
    qc     = ch_qc
    tables = ch_tables
}
