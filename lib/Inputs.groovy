/*
 * Which optional inputs a run actually has, and the steps that follow from them.
 *
 * Nextflow ignores a script-level `params.x = ...` when x is already defined, so a missing input
 * cannot switch its steps off by rewriting the params. Every workflow asks here instead, and the
 * run manifest records these effective values so the reports describe what really ran.
 */
class Inputs {

    // Command-line values arrive as strings, so `--genome null` is the four-character "null".
    static Object unset(v) {
        (v == null || v == 'null' || v.toString().trim() == '') ? null : v
    }

    static boolean hasGenome(Map params) {
        unset(params.genome) || unset(params.genome_accession)
    }

    static boolean hasContaminants(Map params) {
        unset(params.contaminant_fasta) || unset(params.contaminant_dir) || unset(params.contaminant_accessions)
    }

    static boolean longReads(Map params) {
        !params.skip_longread && unset(params.longread_samplesheet)
    }

    // no contaminant genomes = nothing to remove
    static boolean decontam(Map params) {
        !params.skip_decontam && hasContaminants(params)
    }

    // Route A, StringTie (and SQANTI3, which needs both) and the locus gene map all place
    // sequence on the genome
    static boolean genomeGuided(Map params) {
        hasGenome(params) && !params.skip_genome_guided
    }

    static boolean stringtie(Map params) {
        hasGenome(params) && !params.skip_stringtie
    }

    static boolean sqanti(Map params) {
        hasGenome(params) && !params.skip_sqanti
    }

    static boolean locusRebuild(Map params) {
        hasGenome(params) && params.run_locus_rebuild && !params.skip_annotation
    }

    // the headline map cannot be one that was never built
    static String primaryMap(Map params) {
        (params.dge_primary_map == 'locus' && !locusRebuild(params)) ? 'string' : params.dge_primary_map
    }

    // the params as this run used them, for the manifest and the reports
    static Map effective(Map params) {
        [ skip_decontam     : !decontam(params),
          skip_genome_guided: !genomeGuided(params),
          skip_stringtie    : !stringtie(params),
          skip_sqanti       : !sqanti(params),
          run_locus_rebuild : locusRebuild(params),
          dge_primary_map   : primaryMap(params) ]
    }
}
