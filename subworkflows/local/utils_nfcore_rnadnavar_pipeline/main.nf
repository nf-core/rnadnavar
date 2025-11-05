//
// Subworkflow with functionality specific to the nf-core/rnadnavar pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include {SAMPLESHEET_TO_CHANNEL     } from '../samplesheet_to_channel/main.nf'
include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
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
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    versions = Channel.empty()

    // Print version and exit if required and dump pipeline parameters to JSON file
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[1;32m   ___   _  _    _  \033[0;32m  __  _  _    _ \033[1;32m _   _  _    ___  \033[0m
\033[1;32m  | _ \\ | \\| |  /_\\ \033[0;32m |   \\ \\| |  /_\\\033[1;32m\\ \\ / //_\\  | _ \\ \033[0m
\033[1;32m  | _ / | \\` | / _ \\ \033[0;32m| | | \\| | / _ \\\033[1;32m\\ \\ // _ \\ | _ /\033[0m
\033[1;32m  |_|\\_\\|_|\\_|/_/ \\_\\\033[0;32m|___/_|\\_|/_/ \\_\\\033[1;32m\\_//_/ \\_\\|_\\_\\\033[0m
\033[0;35m  nf-core/rnadnavar ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi ? workflow.manifest.doi.tokenize(",").collect { " https://doi.org/${it.trim().replace('https://doi.org/','')}"}.join("\n") : ""}${workflow.manifest.doi ? "\n" : ""} 
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/rnadnavar/blob/master/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // Check input path parameters to see if they exist
    def checkPathParamList = [
        params.input,
        params.fasta,
        params.fasta_fai,
        params.dict,
        params.bwa,
        params.bwamem2,
        params.dragmap,
        params.gtf,
        params.gff,
        params.dbsnp,
        params.dbsnp_tbi,
        params.known_indels,
        params.known_indels_tbi,
        params.germline_resource,
        params.germline_resource_tbi,
        params.intervals,
        params.pon,
        params.pon_tbi,
        params.multiqc_config,
        params.vep_cache,
        params.star_index,
        params.hisat2_index,
        params.whitelist,
    ]


    params.input_restart = retrieveInput((!params.build_only_index && !input), params.step, params.outdir)

    ch_from_samplesheet = params.build_only_index
        ? Channel.empty()
        : input
            ? Channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            : Channel.fromList(samplesheetToList(params.input_restart, "${projectDir}/assets/schema_input.json"))

    SAMPLESHEET_TO_CHANNEL(ch_from_samplesheet)

    emit:
    samplesheet = SAMPLESHEET_TO_CHANNEL.out.input_sample
    versions    = versions
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
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    // Completion email and summary
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
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Check and validate pipeline parameters
def validateInputParameters() {
    genomeExistsError()
}

// Exit pipeline if incorrect --genome key provided
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
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
        "FastQC (Andrews 2010),",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
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
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of the pipeline version used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

// retrieveInput
def retrieveInput(need_input, step, outdir) {
    def input = null
    if (need_input) {
        if (step == 'mapping') {
            error("Can't start ${step} step without samplesheet")
        }
        else if (step == 'markduplicates') {
            log.warn("Using file ${outdir}/csv/mapped.csv")
            input = outdir + "/csv/mapped.csv"
        }
        else if (step == 'prepare_recalibration') {
            log.warn("Using file ${outdir}/csv/markduplicates_no_table.csv")
            input = outdir + "/csv/markduplicates_no_table.csv"
        }
        else if (step == 'recalibrate') {
            log.warn("Using file ${outdir}/csv/markduplicates.csv")
            input = outdir + "/csv/markduplicates.csv"
        }
        else if (step == 'variant_calling') {
            log.warn("Using file ${outdir}/csv/recalibrated.csv")
            input = outdir + "/csv/recalibrated.csv"
        }
        else if (step == 'annotate') {
            log.warn("Using file ${outdir}/csv/variantcalled.csv")
            input = outdir + "/csv/variantcalled.csv"
        }
        else {
            log.warn("Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'")
            error("Unknown step ${step}")
        }
    }
    return input
}
