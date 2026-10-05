#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bioinfo-kaust/hybridmetatranscriptomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/bioinfo-kaust/hybrid-metatranscriptomics
----------------------------------------------------------------------------------------
    Reconstruct a host transcriptome from long + short reads with contaminant reads removed
    BEFORE assembly, then annotate, quantify and test it.

    Written from the Cassiopea andromeda UV experiment and generalised: the species, the
    contaminant set, the experimental design and every tool argument are parameters.
    assets/params_cassiopea_uv.yml reproduces that project; another -params-file runs another.

    The comparison against a published reference deliberately lives OUTSIDE the pipeline —
    see extra/compare_to_published.sh.
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HYBRIDMETATRANSCRIPTOMICS  } from './workflows/hybridmetatranscriptomics'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_hybridmetatranscriptomics_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_hybridmetatranscriptomics_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow BIOINFOKAUST_HYBRIDMETATRANSCRIPTOMICS {

    take:
    samplesheet // channel: samplesheet read in from --input
    assembly    // channel: samplesheet read in from --assembly_input
    longreads   // channel: samplesheet read in from --longread_input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    HYBRIDMETATRANSCRIPTOMICS (
        samplesheet,
        assembly,
        longreads,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report = HYBRIDMETATRANSCRIPTOMICS.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.assembly_input,
        params.longread_input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    BIOINFOKAUST_HYBRIDMETATRANSCRIPTOMICS (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.assembly,
        PIPELINE_INITIALISATION.out.longreads
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        BIOINFOKAUST_HYBRIDMETATRANSCRIPTOMICS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
