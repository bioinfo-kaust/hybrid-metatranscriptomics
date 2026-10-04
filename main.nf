#!/usr/bin/env nextflow
/*
 * hybrid-metatranscriptomics
 * ====================
 * Reconstruct a host transcriptome from long + short reads with contaminant reads removed
 * BEFORE assembly, then annotate, quantify and test it.
 *
 * Written from the Cassiopea andromeda UV experiment (results_host/README.md) and generalised:
 * the species, the contaminant set, the experimental design and every tool argument are
 * parameters. assets/params_cassiopea_uv.yml reproduces that project; another -params-file runs another.
 *
 * The comparison against a published reference deliberately lives OUTSIDE the pipeline —
 * see extra/compare_to_published.sh.
 */

nextflow.enable.dsl = 2

import groovy.json.JsonOutput

include { PREPARE_REFERENCES     } from './subworkflows/local/prepare_references'
include { SHORT_READ_PREPROCESS  } from './subworkflows/local/read_preprocess'
include { SHORT_READ_PREPROCESS as ASSEMBLY_READ_PREPROCESS } from './subworkflows/local/read_preprocess'
include { LONG_READ_PREPROCESS   } from './subworkflows/local/read_preprocess'
include { BUILD_TRANSCRIPTOME    } from './subworkflows/local/build_transcriptome'
include { ANNOTATE               } from './subworkflows/local/annotate'
include { QUANTIFY_AND_DGE       } from './subworkflows/local/quantify_and_dge'
include { BUILD_DECON_INDEX      } from './modules/local/decontam'
include { AGGREGATE_TSV          } from './modules/local/reporting'
include { AGGREGATE_TSV as AGGREGATE_DECONTAM } from './modules/local/reporting'
include { AGGREGATE_TSV as AGGREGATE_ASSEMBLY_READS } from './modules/local/reporting'
include { DECONTAM_MQC           } from './modules/local/reporting'
include { MULTIQC                } from './modules/local/reporting'
include { FIGURES                } from './modules/local/reporting'
include { RUN_SUMMARY            } from './modules/local/reporting'
include { DASHBOARD              } from './modules/local/reporting'

def helpMessage() {
    log.info """
    hybrid-metatranscriptomics ${workflow.manifest.version}

    Usage:
      nextflow run . -profile slurm,conda -params-file my.yml                 # any dataset
      nextflow run . -profile hpc,envmodules -params-file assets/params_cassiopea_uv.yml   # the Cassiopea run (add your read paths)
      nextflow run . -profile hpc,envmodules \\
          --samplesheet short.csv --longread_samplesheet long.csv \\
          --genome_accession GCA_018155075.1 --contaminant_accessions GCA_977109205.1,GCA_965643015.1 \\
          --reference_level Control --contrasts contrasts.csv --busco_lineage cnidaria_odb12 \\
          --eggnog_data_dir /path/to/eggnog/data --slurm_account my_allocation --outdir results

    Core inputs
      --samplesheet            CSV: sample,condition,replicate,fastq_1,fastq_2 (one row per lane)
      --longread_samplesheet   CSV: sample,condition,reads   (null disables the long-read side)
      --assembly_samplesheet   CSV like --samplesheet: extra libraries used only to assemble
      --genome | --genome_accession           reference genome (local FASTA or NCBI accession)
      --contaminant_fasta | --contaminant_dir | --contaminant_accessions
                                              genomes whose reads are removed before assembly
      --outdir                 results directory

    Design
      --condition_column       column in the samplesheet holding the factor  [${params.condition_column}]
      --reference_level        baseline level for the contrasts
      --design                 DESeq2 design formula                          [${params.design}]
      --contrasts              CSV: id,factor,level,reference

    Steps (all default false unless noted)
      --skip_fastp --skip_rrna --skip_decontam --skip_longread
      --skip_genome_guided --skip_reffree --skip_stringtie --skip_denovo
      --skip_annotation --skip_eggnog --skip_interproscan
      --skip_dge --skip_enrichment --skip_dtu (default true) --skip_busco
      --skip_figures --skip_multiqc --skip_dashboard --run_locus_rebuild (default true)

    Run `nextflow run . --help` for this message; see README.md for the full parameter list.
    """.stripIndent()
}

if (params.help) { helpMessage(); exit 0 }

// ---------------------------------------------------------------- samplesheet parsing
def parse_short_samplesheet(path) {
    Channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample) error "samplesheet: every row needs a 'sample' column"
            if (!row[params.condition_column])
                error "samplesheet: missing condition column '${params.condition_column}' in row for ${row.sample}"
            def meta = [ id: row.sample, condition: row[params.condition_column],
                         replicate: row.replicate ?: '' ]
            if (!row.fastq_1 || !row.fastq_2)
                error "samplesheet: ${row.sample} needs both fastq_1 and fastq_2 (this pipeline is paired-end)"
            tuple(meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }
        // several rows per sample = several lanes/units of the same library
        .groupTuple(by: 0)
}

// Assembly-only libraries: same columns, no condition needed; never quantified or tested
def parse_assembly_samplesheet(path) {
    Channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2)
                error "assembly_samplesheet: every row needs sample, fastq_1 and fastq_2"
            def meta = [ id: row.sample, condition: 'assembly_only', replicate: '', assembly_only: true ]
            tuple(meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }
        .groupTuple(by: 0)
}

def parse_long_samplesheet(path) {
    Channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample, condition: row[params.condition_column] ?: 'NA' ]
            tuple(meta, file(row.reads, checkIfExists: true))
        }
}

/*
 * Record the parameters and workflow metadata this run used, so the dashboard can show the
 * arguments each stage really ran with instead of a hard-coded description of the pipeline.
 */
def write_run_manifest() {
    def safe = params.collectEntries { k, v ->
        [(k): (v == null || v instanceof Boolean || v instanceof Number || v instanceof Map || v instanceof List)
              ? v : v.toString()]
    }
    def out = file("${params.outdir}/pipeline_info/run_manifest.json")
    out.parent.mkdirs()
    out.text = JsonOutput.prettyPrint(JsonOutput.toJson([
        workflow: [ run_name    : workflow.runName,
                    start       : workflow.start?.toString(),
                    profile     : workflow.profile,
                    revision    : workflow.revision ?: workflow.manifest.version,
                    nextflow    : workflow.nextflow?.version?.toString(),
                    command_line: workflow.commandLine,
                    outdir      : params.outdir ],
        params: safe
    ]))
    return out
}

workflow {

    log.info """
    ------------------------------------------------------------------
     hybrid-metatranscriptomics ${workflow.manifest.version}
     reference set for quantification : ${params.dge_reference}
     design                           : ${params.design}   (baseline ${params.reference_level ?: 'first level'})
     outdir                           : ${params.outdir}
    ------------------------------------------------------------------
    """.stripIndent()

    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true).first()
    ch_contrasts   = params.contrasts ? Channel.fromPath(params.contrasts, checkIfExists: true).first()
                                      : Channel.fromPath("${projectDir}/assets/NO_FILE_CONTRASTS").first()
    ch_units = parse_short_samplesheet(params.samplesheet)

    // ------------------------------------------------------------------ stage 0
    PREPARE_REFERENCES()

    ch_qc     = Channel.empty()
    ch_tables = Channel.empty()

    // ------------------------------------------------------------------ stage 2 index
    // One combined index: each read's primary alignment is its best placement genome-wide, so a
    // conserved host read that also aligns to a contaminant is still kept.
    ch_decon_index = Channel.empty()
    if (!params.skip_decontam) {
        BUILD_DECON_INDEX(PREPARE_REFERENCES.out.genome, PREPARE_REFERENCES.out.contaminant)
        // both inputs are value channels, so the index is one too and fans out to every sample
        ch_decon_index = BUILD_DECON_INDEX.out.index
    }

    // ------------------------------------------------------------------ stage 1+2 reads
    SHORT_READ_PREPROCESS(ch_units, PREPARE_REFERENCES.out.rrna, ch_decon_index)
    ch_reads  = SHORT_READ_PREPROCESS.out.reads
    ch_qc     = ch_qc.mix(SHORT_READ_PREPROCESS.out.qc)
    ch_tables = ch_tables.mix(SHORT_READ_PREPROCESS.out.tables)

    // Assembly-only libraries: the same preprocessing, then into the two short-read assemblies
    // (StringTie, Trinity) and nowhere else. Their QC tables are published under their own names
    // and kept out of MultiQC and the reports, which describe the experiment's libraries.
    ch_assembly_reads = ch_reads
    ch_assembly_tables = Channel.empty()   // their rRNA/decontam tables, for the dashboard only
    if (params.assembly_samplesheet && params.assembly_samplesheet != 'null') {
        ASSEMBLY_READ_PREPROCESS(parse_assembly_samplesheet(params.assembly_samplesheet),
                                 PREPARE_REFERENCES.out.rrna, ch_decon_index)
        ch_assembly_reads = ch_reads.mix(ASSEMBLY_READ_PREPROCESS.out.reads)
        AGGREGATE_ASSEMBLY_READS(ASSEMBLY_READ_PREPROCESS.out.tables
            .map { name, t -> tuple(name.replace('_per_sample', '') + '_assembly_only', t) })
        ch_assembly_tables = AGGREGATE_ASSEMBLY_READS.out.table.map { name, f -> f }
    }

    ch_long = Channel.empty()
    if (!params.skip_longread && params.longread_samplesheet) {
        LONG_READ_PREPROCESS(parse_long_samplesheet(params.longread_samplesheet),
                             PREPARE_REFERENCES.out.rrna,
                             PREPARE_REFERENCES.out.genome,
                             PREPARE_REFERENCES.out.contaminant)
        ch_long   = LONG_READ_PREPROCESS.out.seqs
        ch_qc     = ch_qc.mix(LONG_READ_PREPROCESS.out.qc)
        ch_tables = ch_tables.mix(LONG_READ_PREPROCESS.out.tables)
    }

    // ------------------------------------------------------------------ stage 3 reference
    BUILD_TRANSCRIPTOME(ch_long, ch_assembly_reads,
                        PREPARE_REFERENCES.out.genome,
                        PREPARE_REFERENCES.out.contaminant,
                        PREPARE_REFERENCES.out.star_index)
    ch_qc = ch_qc.mix(BUILD_TRANSCRIPTOME.out.qc)

    // ------------------------------------------------------------------ stage 4 annotation
    ch_annotation = Channel.empty()
    if (!params.skip_annotation) {
        ANNOTATE(BUILD_TRANSCRIPTOME.out.selected,
                 BUILD_TRANSCRIPTOME.out.selected_tx2gene,
                 PREPARE_REFERENCES.out.diamond_dbs)
        ch_annotation = ANNOTATE.out.table
    }

    // ------------------------------------------------------------------ stage 5+6
    QUANTIFY_AND_DGE(ch_reads,
                     BUILD_TRANSCRIPTOME.out.selected,
                     BUILD_TRANSCRIPTOME.out.selected_tx2gene,
                     ch_annotation,
                     PREPARE_REFERENCES.out.genome,
                     BUILD_TRANSCRIPTOME.out.collapsed_gtf,
                     ch_samplesheet,
                     ch_contrasts)

    // Salmon's own logs, so MultiQC has a quantification section: mapping rate and library type
    // per sample lived only in the bespoke dashboard before this.
    ch_qc = ch_qc.mix(QUANTIFY_AND_DGE.out.quants.flatten())
    // BUSCO writes a short_summary MultiQC can parse; only the one-line grep reached reporting
    if (!params.skip_busco) {
        // only the labelled top-level summaries: BUSCO also writes run_<lineage>/short_summary.txt
        // and .json, whose names are identical across rungs and collide when staged together
        ch_qc = ch_qc.mix(BUILD_TRANSCRIPTOME.out.busco_dir
                            .flatten().filter { f -> f.name ==~ /short_summary\.specific\..*\.txt/ })
    }

    // ------------------------------------------------------------------ reporting
    AGGREGATE_TSV(ch_tables)
    ch_qc_tables = AGGREGATE_TSV.out.table.map { name, f -> f }.collect()

    if (!params.skip_decontam) {
        DECONTAM_MQC(AGGREGATE_TSV.out.table.filter { n, f -> n == 'decontam_per_sample' }.map { n, f -> f })
        ch_qc = ch_qc.mix(DECONTAM_MQC.out.mqc)
    }

    if (!params.skip_multiqc) {
        def mqc_cfg = params.multiqc_config ? file(params.multiqc_config)
                                            : file("${projectDir}/assets/NO_FILE_MQC")
        MULTIQC(ch_qc.collect(), mqc_cfg)
    }

    ch_figures = Channel.empty()
    if (!params.skip_figures && !params.skip_dge) {
        // One figure set per gene map, not only the string-stripped one: the diagnostics (MA,
        // dispersion, p-value shape) and the volcanoes describe a particular map's DESeq2 fit, and
        // the primary map is whichever dge_primary_map names. Each map's DGE dir is joined to its
        // own enrichment dir by label; a map with no enrichment gets the placeholder.
        ch_fig_in = QUANTIFY_AND_DGE.out.dge
            .mix(QUANTIFY_AND_DGE.out.dge_corset, QUANTIFY_AND_DGE.out.dge_locus)
            .join(QUANTIFY_AND_DGE.out.enrichment
                    .mix(QUANTIFY_AND_DGE.out.enrich_corset, QUANTIFY_AND_DGE.out.enrich_locus),
                  remainder: true)
            .map { label, dge_dir, enr_dir -> tuple(label, dge_dir, enr_dir ?: file("${projectDir}/assets/NO_FILE_ENR")) }
        // distinct placeholders: two inputs staged under the same name would collide
        FIGURES(ch_fig_in,
                ch_annotation.ifEmpty { tuple('none', file("${projectDir}/assets/NO_FILE_ANN")) }.first(),
                ch_qc_tables,
                BUILD_TRANSCRIPTOME.out.busco.map { l, f -> f }.collect().ifEmpty([]),
                // value channels from here on, so every map's FIGURES task receives them
                BUILD_TRANSCRIPTOME.out.ladder.first())
        // the dashboard embeds the diagnostic panels, so it needs the figure set staged
        ch_figures = FIGURES.out.figures
    }

    // one place that carries every headline number with the file it came from
    ch_ann_summary = params.skip_annotation ? Channel.empty() : ANNOTATE.out.summary

    // everything the reports are built from, in one place so both consumers agree
    ch_report_inputs = ch_qc_tables
        .mix(BUILD_TRANSCRIPTOME.out.ladder)
        .mix(ch_ann_summary)
        .mix(QUANTIFY_AND_DGE.out.mapping)
        // the DIRECTORIES, not the files inside them: every map writes a DE_summary.csv and a
        // filtering_summary.csv, and three files of the same name cannot stage side by side
        .mix(QUANTIFY_AND_DGE.out.dge.map { l, d -> d }.ifEmpty([]))
        .mix(QUANTIFY_AND_DGE.out.dge_locus.map { l, d -> d }.ifEmpty([]))
        .mix(QUANTIFY_AND_DGE.out.locus_report.ifEmpty([]))
        .mix(QUANTIFY_AND_DGE.out.expr_filter.ifEmpty([]))
        .mix(QUANTIFY_AND_DGE.out.expr_qc.ifEmpty([]))
        .mix(QUANTIFY_AND_DGE.out.corset.ifEmpty([]))
        .mix(QUANTIFY_AND_DGE.out.dge_corset.map { l, d -> d }.ifEmpty([]))
        .mix(BUILD_TRANSCRIPTOME.out.sqanti.ifEmpty([]))

    if (!params.skip_dashboard) {
        // ch_report_inputs already carries the DGE directories and the run's tables; only the
        // extras the dashboard needs on top are mixed here. Staging the same directory twice is
        // an input-name collision, not a no-op, and it fails the task.
        DASHBOARD(
            ch_report_inputs
                .mix(BUILD_TRANSCRIPTOME.out.busco.map { l, f -> f }.ifEmpty([]))
                .mix(params.skip_decontam ? Channel.empty() : BUILD_DECON_INDEX.out.log)
                .mix(ch_annotation.map { l, f -> f }.ifEmpty([]))
                .mix(QUANTIFY_AND_DGE.out.enrichment.map { l, d -> d }.ifEmpty([]))
                .mix(QUANTIFY_AND_DGE.out.enrich_corset.map { l, d -> d }.ifEmpty([]))
                // all three, not just two: the dashboard draws one enrichment box per gene map,
                // and the primary map is whichever dge_primary_map names. Staging only the string
                // and corset arms left the locus arm -- the primary one in this project -- out of
                // the report entirely, so the Enrichment tab described a different map from the
                // one the headline DE numbers came from.
                .mix(QUANTIFY_AND_DGE.out.enrich_locus.map { l, d -> d }.ifEmpty([]))
                .mix(ch_figures.ifEmpty([]))
                .mix(ch_assembly_tables)
                .mix(QUANTIFY_AND_DGE.out.gene_map_files.ifEmpty([]))
                .unique()
                .collect(),
            Channel.fromPath(write_run_manifest()))
    }

    RUN_SUMMARY(ch_report_inputs.unique().collect())

    // Provenance: every process writes a versions.yml and publishes it to the `versions` topic,
    // so this collects the whole run's tool versions without wiring 56 output channels by hand.
    // Sorted and deduplicated, because a process that runs per sample emits one file per task.
    Channel.topic('versions')
        .map { f -> f.text }
        .unique()
        .collectFile(name: 'software_versions.yml', sort: true,
                     storeDir: "${params.outdir}/pipeline_info")
}

workflow.onComplete {
    log.info(workflow.success
        ? "\nDone. Results: ${params.outdir}\nStart with dashboard/index.html (or run_summary.txt for the plain-text version).\n"
        : "\nFailed: ${workflow.errorMessage}\nResume with -resume once the cause is fixed.\n")
}
