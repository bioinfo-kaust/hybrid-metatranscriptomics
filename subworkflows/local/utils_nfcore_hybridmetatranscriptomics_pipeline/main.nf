//
// Subworkflow with functionality specific to the bioinfo-kaust/hybridmetatranscriptomics pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    assembly_input    //  string: Path to the assembly-only samplesheet (or null)
    longread_input    //  string: Path to the long-read samplesheet (or null)
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --longread_input long.csv --genome genome.fa --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Optional inputs: say up front which steps a missing input switches off
    //
    checkOptionalInputs()

    //
    // Create channels from the samplesheets
    //
    // Several rows per sample = several lanes of one library, merged before anything else.
    ch_samplesheet = shortReadChannel(input, false)
    ch_assembly    = isUnset(assembly_input) ? channel.empty() : shortReadChannel(assembly_input, true)
    ch_longreads   = runLongReads()
        ? channel.fromList(longReadList(longread_input))
        : channel.empty()

    emit:
    samplesheet = ch_samplesheet   // [ meta, [ r1_L1, r2_L1, r1_L2, r2_L2, ... ] ]
    assembly    = ch_assembly      // same shape; assembly-only libraries
    longreads   = ch_longreads     // [ meta, reads ]
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        resultsPointer(outdir)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Which optional inputs a run actually has, and the steps that follow from them.
//
// A missing input cannot switch its steps off by rewriting params (a script-level `params.x = ...`
// is ignored when x is already defined), so every workflow asks these instead, and the run
// manifest records the effective values so the reports describe what really ran.
//

// Command-line values arrive as strings, so `--genome null` is the four-character "null"
def isUnset(v) {
    v == null || v == 'null' || v.toString().trim() == ''
}

def hasGenome() {
    !isUnset(params.genome) || !isUnset(params.genome_accession)
}

def hasContaminants() {
    !isUnset(params.contaminant_fasta) || !isUnset(params.contaminant_dir) || !isUnset(params.contaminant_accessions)
}

def runLongReads() {
    !params.skip_longread && !isUnset(params.longread_input)
}

// no contaminant genomes = nothing to remove
def runDecontam() {
    !params.skip_decontam && hasContaminants()
}

// Route A, StringTie (and SQANTI3, which needs both) and the locus gene map all place sequence
// on the genome
def runGenomeGuided() {
    hasGenome() && !params.skip_genome_guided
}

def runStringtie() {
    hasGenome() && !params.skip_stringtie
}

def runSqanti() {
    hasGenome() && !params.skip_sqanti
}

def runLocusRebuild() {
    hasGenome() && params.run_locus_rebuild && !params.skip_annotation
}

// the headline map cannot be one that was never built
def primaryGeneMap() {
    (params.dge_primary_map == 'locus' && !runLocusRebuild()) ? 'string' : params.dge_primary_map
}

// the params as this run used them, for the manifest and the reports
def effectiveParams() {
    [ skip_decontam     : !runDecontam(),
      skip_genome_guided: !runGenomeGuided(),
      skip_stringtie    : !runStringtie(),
      skip_sqanti       : !runSqanti(),
      run_locus_rebuild : runLocusRebuild(),
      dge_primary_map   : primaryGeneMap() ]
}

//
// Startup checks for the optional inputs and the one tool without a container-free provider
//
def checkOptionalInputs() {
    if (!hasContaminants() && !params.skip_decontam) {
        log.warn "No contaminant genomes given: contaminant removal is skipped."
    }
    if (!hasGenome()) {
        log.warn "No genome given: skipping the genome-guided long-read route, StringTie, SQANTI3 and the " +
                 "locus gene map${params.dge_primary_map == 'locus' ? " (primary map falls back to 'string')" : ''}; " +
                 "salmon runs without decoys." +
                 (runDecontam() ? " Decontamination removes every read that places on a contaminant." : '')
    }
    def has_long = runLongReads()
    if (!(has_long && runGenomeGuided()) && !runStringtie() && params.skip_denovo && !(has_long && !params.skip_reffree)) {
        error "No assembly route left to build a transcriptome from: give a --genome, a long-read " +
              "samplesheet, or leave Trinity on (--skip_denovo false)."
    }
    if (!params.skip_eggnog && !params.skip_annotation && isUnset(params.eggnog_data_dir)) {
        error "eggNOG-mapper needs --eggnog_data_dir (or run with --skip_eggnog)."
    }

    // corset has no environment module: outside a container or conda it must come from
    // --corset_bin, --extra_path or the PATH. Look now, not when CORSET fails hours into the run.
    def provisioned = workflow.containerEngine || workflow.profile.tokenize(',').intersect(['conda', 'mamba'])
    if (!params.skip_corset && !provisioned && !workflow.stubRun) {
        def dirs = [params.corset_bin, params.extra_path].findAll { d -> !isUnset(d) } +
                   (System.getenv('PATH') ?: '').tokenize(File.pathSeparator)
        if (!isUnset(params.corset_bin) && !file("${params.corset_bin}/corset").exists()) {
            error "--corset_bin ${params.corset_bin} has no 'corset' executable."
        }
        if (!dirs.any { d -> file("${d}/corset").exists() }) {
            error "corset not found (no environment module provides it): set --corset_bin to the " +
                  "directory holding it, use a container or conda profile, or run with --skip_corset."
        }
    }
}

//
// Short-read samplesheets: --input (the experiment) and --assembly_input (assembly-only libraries)
//
def shortReadChannel(sheet, assembly_only) {
    // the factor column is a parameter, so it is read from the sheet rather than the schema
    def rows = file(sheet, checkIfExists: true).splitCsv(header: true)
    def condition = rows.collectEntries { row -> [(row.sample): row[params.condition_column]] }
    if (!assembly_only) {
        def missing = condition.findAll { _id, c -> !c }.keySet()
        if (missing) {
            error "samplesheet: missing condition column '${params.condition_column}' for ${missing.join(', ')}"
        }
    }
    return channel
        .fromList(samplesheetToList(sheet, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2 ->
            def m = assembly_only
                ? [ id: meta.id, condition: 'assembly_only', replicate: '', assembly_only: true, single_end: false ]
                : [ id: meta.id, condition: condition[meta.id], replicate: meta.replicate ?: '', single_end: false ]
            [ m.id, m, [ fastq_1, fastq_2 ] ]
        }
        .groupTuple()
        .map { _id, metas, lanes ->
            // lanes in file-name order, R1/R2 interleaved: the order `cat` concatenates them in
            [ metas[0], lanes.sort { a, b -> a[0].name <=> b[0].name }.flatten() ]
        }
}

//
// Long-read samplesheet: sample,condition,reads
//
def longReadList(sheet) {
    def rows = file(sheet, checkIfExists: true).splitCsv(header: true)
    def condition = rows.collectEntries { row -> [(row.sample): row[params.condition_column]] }
    return samplesheetToList(sheet, "${projectDir}/assets/schema_longread.json").collect { meta, reads ->
        [ [ id: meta.id, condition: condition[meta.id] ?: 'NA', single_end: true ], reads ]
    }
}

//
// Record the parameters and workflow metadata this run used, so the dashboard can show the
// arguments each stage really ran with instead of a hard-coded description of the pipeline
//
def writeRunManifest(outdir) {
    // the effective values, not the requested ones: a step switched off by a missing input
    // must read as off in the reports
    def safe = (params + effectiveParams()).collectEntries { k, v ->
        [(k): (v == null || v instanceof Boolean || v instanceof Number || v instanceof Map || v instanceof List)
              ? v : v.toString()]
    }
    def out = file("${outdir}/pipeline_info/run_manifest.json")
    out.parent.mkdirs()
    out.text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        workflow: [ run_name    : workflow.runName,
                    start       : workflow.start?.toString(),
                    profile     : workflow.profile,
                    revision    : workflow.revision ?: workflow.manifest.version,
                    nextflow    : workflow.nextflow?.version?.toString(),
                    command_line: workflow.commandLine,
                    outdir      : outdir ],
        params: safe
    ]))
    return out
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Where to start reading the results (inside onComplete `workflow` is not the run metadata)
//
def resultsPointer(outdir) {
    if (workflow.success) {
        log.info "Results: ${outdir}\nStart with dashboard/index.html (or run_summary.txt for the plain-text version)."
    }
}
