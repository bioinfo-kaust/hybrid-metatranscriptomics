/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPARE_REFERENCES                                } from '../subworkflows/local/prepare_references'
include { SHORT_READ_PREPROCESS                             } from '../subworkflows/local/short_read_preprocess'
include { SHORT_READ_PREPROCESS as ASSEMBLY_READ_PREPROCESS } from '../subworkflows/local/short_read_preprocess'
include { LONG_READ_PREPROCESS                              } from '../subworkflows/local/long_read_preprocess'
include { BUILD_TRANSCRIPTOME                               } from '../subworkflows/local/build_transcriptome'
include { EXPRESSION_FILTER                                 } from '../subworkflows/local/expression_filter'
include { ANNOTATE                                          } from '../subworkflows/local/annotate'
include { QUANTIFY_AND_DGE                                  } from '../subworkflows/local/quantify_and_dge'
include { BUILD_DECON_INDEX                                 } from '../modules/local/build_decon_index'
include { AGGREGATE_TSV                                     } from '../modules/local/aggregate_tsv'
include { AGGREGATE_TSV as AGGREGATE_ASSEMBLY_READS         } from '../modules/local/aggregate_tsv'
include { DECONTAM_MQC                                      } from '../modules/local/decontam_mqc'
include { FIGURES                                           } from '../modules/local/figures'
include { RUN_SUMMARY                                       } from '../modules/local/run_summary'
include { DASHBOARD                                         } from '../modules/local/dashboard'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                            } from '../subworkflows/local/utils_nfcore_hybridmetatranscriptomics_pipeline'
include { runDecontam; runLongReads; primaryGeneMap; writeRunManifest; isUnset } from '../subworkflows/local/utils_nfcore_hybridmetatranscriptomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HYBRIDMETATRANSCRIPTOMICS {

    take:
    ch_samplesheet   // channel: [ meta, [ reads ] ] read in from --input
    ch_assembly      // channel: [ meta, [ reads ] ] read in from --assembly_input (or empty)
    ch_longreads     // channel: [ meta, reads ] read in from --longread_input (or empty)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions      = channel.empty()
    def ch_multiqc_files = channel.empty()

    def expr_rule = params.expr_filter_tpm > 0
        ? "TPM >= ${params.expr_filter_tpm} in >= ${params.expr_filter_samples} libraries"
        : 'off (the whole set is annotated and tested)'
    log.info """
    ------------------------------------------------------------------
     assembled reference set : ${params.dge_reference}
     expression filter       : ${expr_rule}
                               annotation, DGE and enrichment all run on what it keeps
     design                  : ${params.design}   (baseline ${params.reference_level ?: 'first level'})
     outdir                  : ${outdir}
    ------------------------------------------------------------------
    """.stripIndent()

    ch_sheet     = channel.value(file(params.input, checkIfExists: true))
    ch_contrasts = channel.value(isUnset(params.contrasts)
        ? file("${projectDir}/assets/NO_FILE_CONTRASTS", checkIfExists: true)
        : file(params.contrasts, checkIfExists: true))

    // ------------------------------------------------------------------ stage 0
    PREPARE_REFERENCES()

    ch_tables = channel.empty()

    // ------------------------------------------------------------------ stage 2 index
    // One combined index: each read's primary alignment is its best placement genome-wide, so a
    // conserved host read that also aligns to a contaminant is still kept.
    ch_decon_index = channel.empty()
    ch_decon_log   = channel.empty()
    if (runDecontam()) {
        BUILD_DECON_INDEX(PREPARE_REFERENCES.out.genome, PREPARE_REFERENCES.out.contaminant, params.contaminant_prefix)
        // both inputs are value channels, so the index is one too and fans out to every sample
        ch_decon_index = BUILD_DECON_INDEX.out.index
        ch_decon_log   = BUILD_DECON_INDEX.out.log
    }

    // ------------------------------------------------------------------ stage 1+2 reads
    SHORT_READ_PREPROCESS(ch_samplesheet, PREPARE_REFERENCES.out.rrna, ch_decon_index)
    ch_reads         = SHORT_READ_PREPROCESS.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(SHORT_READ_PREPROCESS.out.qc)
    ch_tables        = ch_tables.mix(SHORT_READ_PREPROCESS.out.tables)

    // Assembly-only libraries: the same preprocessing, then into the two short-read assemblies
    // (StringTie, Trinity) and nowhere else. Their QC tables are published under their own names
    // and kept out of MultiQC and the reports, which describe the experiment's libraries.
    ASSEMBLY_READ_PREPROCESS(ch_assembly, PREPARE_REFERENCES.out.rrna, ch_decon_index)
    ch_assembly_reads = ch_reads.mix(ASSEMBLY_READ_PREPROCESS.out.reads)
    AGGREGATE_ASSEMBLY_READS(ASSEMBLY_READ_PREPROCESS.out.tables
        .map { name, t -> [ name.replace('_per_sample', '') + '_assembly_only', t ] })
    // their rRNA/decontam tables, for the dashboard only
    ch_assembly_tables = AGGREGATE_ASSEMBLY_READS.out.table.map { _name, f -> f }

    ch_long = channel.empty()
    if (runLongReads()) {
        LONG_READ_PREPROCESS(ch_longreads,
                             PREPARE_REFERENCES.out.rrna,
                             PREPARE_REFERENCES.out.genome,
                             PREPARE_REFERENCES.out.contaminant)
        ch_long          = LONG_READ_PREPROCESS.out.seqs
        ch_multiqc_files = ch_multiqc_files.mix(LONG_READ_PREPROCESS.out.qc)
        ch_tables        = ch_tables.mix(LONG_READ_PREPROCESS.out.tables)
    }

    // ------------------------------------------------------------------ stage 3 reference
    BUILD_TRANSCRIPTOME(ch_long, ch_assembly_reads,
                        PREPARE_REFERENCES.out.genome,
                        PREPARE_REFERENCES.out.contaminant,
                        PREPARE_REFERENCES.out.star_index)
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_TRANSCRIPTOME.out.qc)

    // ------------------------------------------------------------------ stage 3b expression filter
    /*
     * The analysed set is decided here, before anything is annotated or tested. The ladder tip is
     * a cross-source union with no expression evidence applied at any rung, so one salmon pass
     * measures what the experiment actually saw and the reference is cut to the survivors; every
     * step after this — annotation, quantification, DE, enrichment, figures, reports — describes
     * that same set, so no number has to be qualified with which reference it came from.
     *
     * The ASSEMBLY is untouched: BUILD_TRANSCRIPTOME publishes every rung with its own seqkit and
     * BUSCO numbers, so the completeness quoted for the transcriptome is still a property of the
     * full deduplicated set. EXPRESSION_FILTER adds the same two measurements for the filtered set.
     *
     * expr_filter_tpm = 0 switches the filter off: the ladder tip itself is annotated and tested,
     * and the run has one salmon pass either way.
     */
    ch_ref_fasta   = BUILD_TRANSCRIPTOME.out.selected
    ch_ref_tx2gene = BUILD_TRANSCRIPTOME.out.selected_tx2gene
    ch_expr_log    = channel.empty()
    ch_expr_qc     = channel.empty()
    // no filter -> every gene in the reference is DE-eligible, signalled by the sentinel
    ch_expr_tx     = channel.value(file("${projectDir}/assets/NO_FILE_EXPRTX", checkIfExists: true))
    if (params.expr_filter_tpm > 0) {
        EXPRESSION_FILTER(ch_reads, ch_ref_fasta, ch_ref_tx2gene, PREPARE_REFERENCES.out.genome)
        ch_ref_fasta   = EXPRESSION_FILTER.out.fasta
        ch_ref_tx2gene = EXPRESSION_FILTER.out.tx2gene
        ch_expr_log    = EXPRESSION_FILTER.out.filter_log
        ch_expr_qc     = EXPRESSION_FILTER.out.qc
        // the reference IS the expressed set, so this restriction is an identity; it stays wired
        // so that a second way into the reference cannot silently widen what gets tested
        ch_expr_tx     = EXPRESSION_FILTER.out.expressed_tx.first()
    }

    // ------------------------------------------------------------------ stage 4 annotation
    ch_annotation  = channel.empty()
    ch_ann_summary = channel.empty()
    // no annotation -> the gene-level FASTA is nucleotide only, signalled by the sentinel
    ch_pep         = channel.value(file("${projectDir}/assets/NO_FILE_PEP", checkIfExists: true))
    if (!params.skip_annotation) {
        ANNOTATE(ch_ref_fasta,
                 ch_ref_tx2gene,
                 PREPARE_REFERENCES.out.diamond_dbs)
        ch_annotation  = ANNOTATE.out.table
        ch_ann_summary = ANNOTATE.out.summary
        ch_pep         = ANNOTATE.out.pep.map { _meta, pep -> pep }.first()
    }

    // ------------------------------------------------------------------ stage 5+6
    QUANTIFY_AND_DGE(ch_reads,
                     ch_ref_fasta,
                     ch_ref_tx2gene,
                     ch_annotation,
                     ch_pep,
                     ch_expr_tx,
                     PREPARE_REFERENCES.out.genome,
                     BUILD_TRANSCRIPTOME.out.collapsed_gtf,
                     ch_sheet,
                     ch_contrasts)
    ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_AND_DGE.out.qc)

    // ------------------------------------------------------------------ reporting
    AGGREGATE_TSV(ch_tables)
    ch_qc_tables = AGGREGATE_TSV.out.table.map { _name, f -> f }.collect()

    if (runDecontam()) {
        DECONTAM_MQC(AGGREGATE_TSV.out.table.filter { n, _f -> n == 'decontam_per_sample' }.map { _n, f -> f })
        ch_multiqc_files = ch_multiqc_files.mix(DECONTAM_MQC.out.mqc)
    }

    ch_figures = channel.empty()
    if (!params.skip_figures && !params.skip_dge) {
        // One figure set per gene map, not only the string-stripped one: the diagnostics (MA,
        // dispersion, p-value shape) and the volcanoes describe a particular map's DESeq2 fit, and
        // the primary map is whichever dge_primary_map names. Each map's DGE dir is joined to its
        // own enrichment dir; a map with no enrichment gets the placeholder.
        ch_fig_in = QUANTIFY_AND_DGE.out.dge
            .mix(QUANTIFY_AND_DGE.out.dge_corset, QUANTIFY_AND_DGE.out.dge_locus)
            .join(QUANTIFY_AND_DGE.out.enrichment
                    .mix(QUANTIFY_AND_DGE.out.enrich_corset, QUANTIFY_AND_DGE.out.enrich_locus),
                  remainder: true)
            .map { meta, dge_dir, enr_dir -> [ meta, dge_dir, enr_dir ?: file("${projectDir}/assets/NO_FILE_ENR") ] }
        // distinct placeholders: two inputs staged under the same name would collide
        FIGURES(ch_fig_in,
                ch_annotation.ifEmpty([ [id: 'none'], file("${projectDir}/assets/NO_FILE_ANN") ]).first(),
                ch_qc_tables,
                BUILD_TRANSCRIPTOME.out.busco.map { _meta, f -> f }.collect().ifEmpty([]),
                // value channels from here on, so every map's FIGURES task receives them
                BUILD_TRANSCRIPTOME.out.ladder.first())
        // the dashboard embeds the diagnostic panels, so it needs the figure set staged
        ch_figures = FIGURES.out.figures
    }

    // everything the reports are built from, in one place so both consumers agree
    ch_report_inputs = ch_qc_tables
        .mix(BUILD_TRANSCRIPTOME.out.ladder)
        .mix(ch_ann_summary)
        .mix(QUANTIFY_AND_DGE.out.mapping)
        // the DIRECTORIES, not the files inside them: every map writes a DE_summary.csv and a
        // filtering_summary.csv, and three files of the same name cannot stage side by side
        .mix(QUANTIFY_AND_DGE.out.dge.map { _meta, d -> d })
        .mix(QUANTIFY_AND_DGE.out.dge_locus.map { _meta, d -> d })
        .mix(QUANTIFY_AND_DGE.out.locus_report)
        .mix(ch_expr_log)
        .mix(ch_expr_qc)
        .mix(QUANTIFY_AND_DGE.out.corset)
        .mix(QUANTIFY_AND_DGE.out.dge_corset.map { _meta, d -> d })
        .mix(BUILD_TRANSCRIPTOME.out.sqanti)

    if (!params.skip_dashboard) {
        // ch_report_inputs already carries the DGE directories and the run's tables; only the
        // extras the dashboard needs on top are mixed here. Staging the same directory twice is
        // an input-name collision, not a no-op, and it fails the task.
        DASHBOARD(
            ch_report_inputs
                .mix(BUILD_TRANSCRIPTOME.out.busco.map { _meta, f -> f })
                .mix(ch_decon_log)
                .mix(ch_annotation.map { _meta, f -> f })
                .mix(QUANTIFY_AND_DGE.out.enrichment.map { _meta, d -> d })
                .mix(QUANTIFY_AND_DGE.out.enrich_corset.map { _meta, d -> d })
                // all three: the dashboard draws one enrichment box per gene map, and the primary
                // map is whichever dge_primary_map names
                .mix(QUANTIFY_AND_DGE.out.enrich_locus.map { _meta, d -> d })
                .mix(ch_figures)
                .mix(ch_assembly_tables)
                .mix(QUANTIFY_AND_DGE.out.gene_map_files)
                .mix(QUANTIFY_AND_DGE.out.gene_rep)
                .unique()
                .collect(),
            channel.value(writeRunManifest(outdir)),
            params.dashboard_title ?: "Transcriptome run — ${params.dge_reference}")
    }

    RUN_SUMMARY(ch_report_inputs.unique().collect(), primaryGeneMap())

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'hybridmetatranscriptomics_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'hybridmetatranscriptomics'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        ch_multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }
    }

    emit:
    multiqc_report = ch_multiqc_report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
