/*
 * Stage 1+2 (long reads) — normalise -> rRNA depletion -> pool -> contaminant classification.
 */

include { LONGREAD_PREP       } from '../../../modules/local/longread_prep'
include { LONGREAD_MEMBERSHIP } from '../../../modules/local/longread_membership'
include { BBDUK_RRNA_SE       } from '../../../modules/local/bbduk_rrna_se'
include { LONGREAD_POOL       } from '../../../modules/local/longread_pool'
include { DECONTAM_LONG       } from '../../../modules/local/decontam_long'
include { isUnset; runDecontam } from '../utils_nfcore_hybridmetatranscriptomics_pipeline'

workflow LONG_READ_PREPROCESS {

    take:
    ch_long           // [ meta, reads ]
    ch_rrna           // path (value)
    ch_genome         // path (value)
    ch_contaminant    // path (value)

    main:
    ch_qc     = channel.empty()
    ch_tables = channel.empty()

    LONGREAD_PREP(ch_long)
    ch_seqs = LONGREAD_PREP.out.seqs

    // which sequence came from which sample, taken before anything is filtered so the per-sample
    // input counts are the real ones. An external map (--longread_membership) covers the case
    // where the long reads were already pooled before the pipeline saw them.
    ch_membership = !isUnset(params.longread_membership)
        ? channel.fromPath(params.longread_membership, checkIfExists: true).collect()
        : LONGREAD_MEMBERSHIP(ch_seqs).tsv.collect().ifEmpty([ file("${projectDir}/assets/NO_FILE_MEMBERSHIP") ])

    ch_rrna_hits = channel.fromPath("${projectDir}/assets/NO_FILE_RRNA").collect()
    if (!params.skip_rrna) {
        BBDUK_RRNA_SE(ch_seqs, ch_rrna)
        ch_seqs      = BBDUK_RRNA_SE.out.seqs
        ch_qc        = ch_qc.mix(BBDUK_RRNA_SE.out.stats)
        ch_tables    = ch_tables.mix(BBDUK_RRNA_SE.out.summary.collect().map { t -> [ 'longread_rrna', t ] })
        ch_rrna_hits = BBDUK_RRNA_SE.out.removed.collect()
    }

    // sample -> condition, for the per-sample table
    ch_conditions = ch_long
        .map { meta, _reads -> "${meta.id}\t${meta.condition}\n" }
        .collectFile(name: 'longread_conditions.tsv', sort: true)

    // The reference is built from all long reads together: these libraries are usually few and
    // unreplicated, and isoform discovery benefits from the pooled depth.
    ch_pooled = LONGREAD_POOL(ch_seqs.map { _meta, seqs -> seqs }.collect()).seqs

    if (runDecontam()) {
        DECONTAM_LONG(ch_pooled, ch_genome, ch_contaminant, ch_membership, ch_rrna_hits, ch_conditions)
        ch_pooled = DECONTAM_LONG.out.seqs
        ch_tables = ch_tables
            .mix(DECONTAM_LONG.out.summary.map { t -> [ 'longread_decontam', [t] ] })
            .mix(DECONTAM_LONG.out.per_sample.map { t -> [ 'longread_per_sample', [t] ] })
    }

    emit:
    seqs   = ch_pooled   // [ [id: 'longreads'], fasta ]
    qc     = ch_qc
    tables = ch_tables
}
