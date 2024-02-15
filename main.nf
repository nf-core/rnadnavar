#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnadnavar
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnadnavar
    Website: https://nf-co.re/rnadnavar
    Slack  : https://nfcore.slack.com/channels/rnadnavar
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.bwa                  = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.bwamem2              = WorkflowMain.getGenomeAttribute(params, 'bwamem2')
params.fasta                = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai            = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.dict                 = WorkflowMain.getGenomeAttribute(params, 'dict')
params.gtf                  = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.gff                  = WorkflowMain.getGenomeAttribute(params, 'gff')
params.exon_bed             = WorkflowMain.getGenomeAttribute(params, 'exon_bed')
params.star_index           = WorkflowMain.getGenomeAttribute(params, 'star')
params.dbsnp                = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_tbi            = WorkflowMain.getGenomeAttribute(params, 'dbsnp_tbi')
params.known_indels         = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_tbi     = WorkflowMain.getGenomeAttribute(params, 'known_indels_tbi')
params.vep_cache_version    = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.vep_genome           = WorkflowMain.getGenomeAttribute(params, 'vep_genome')
params.vep_species          = WorkflowMain.getGenomeAttribute(params, 'vep_species')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ALTERNATIVE INPUT FILE ON RESTART
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.input_restart = WorkflowRnadnavar.retrieveInput(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def String command = "nextflow run ${workflow.manifest.name} --input samplesheet.csv --genome GATK.GRCh38 -profile docker --outdir results"
    log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}

WorkflowMain.initialise(workflow, params, log, args)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNADNAVAR } from './workflows/rnadnavar'

//
// WORKFLOW: Run main nf-core/rnadnavar analysis pipeline
//
workflow NFCORE_RNADNAVAR {
    RNADNAVAR ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
workflow {
    NFCORE_RNADNAVAR ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
