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
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
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
    // Build text sections for different pipeline components
    def text_preprocessing = [
        "Raw read quality control was performed with FastQC (Andrews 2010)",
        params.trimmer == "fastp" ? "and preprocessing with fastp (Chen et al. 2018)." : "."
    ].join(' ').trim()

    def text_alignment = [
        "Read alignment was performed with",
        params.aligner == "star" ? "STAR (Dobin et al. 2013)" : "",
        params.aligner == "bwamem" ? "BWA-MEM (Li 2013)" : "",
        params.aligner == "bwamem2" ? "BWA-MEM2 (Vasimuddin et al. 2019)" : "",
        params.aligner == "dragmap" ? "DragMap" : "",
        params.tools.contains("realignment") ? "HISAT2 (Kim et al. 2019)" : "",
        "for both DNA and RNA samples."
    ].join(' ').trim()

    def text_alignment_processing = [
        "Alignment files were processed with SAMtools (Li et al. 2009)",
        params.run_mosdepth ? "and coverage analysis with Mosdepth (Pedersen & Quinlan 2018)." : "."
    ].join(' ').trim()

    def text_gatk_preprocessing = [
        params.run_gatk_preprocessing ? "GATK preprocessing including base quality score recalibration and duplicate marking was performed (McKenna et al. 2010)." : ""
    ].join(' ').trim()

    def text_variant_calling = [
        "Variant calling was performed using",
        params.tools.contains("mutect2") ? "GATK Mutect2 (McKenna et al. 2010)," : "",
        params.tools.contains("strelka2") ? "Strelka2 (Kim et al. 2018)," : "",
        params.tools.contains("sage") ? "SAGE," : "",
        "with consensus calling when multiple callers were used."
    ].join(' ').trim().replaceAll(",\\s*with", " with").replaceAll(",\$", ".")

    def text_variant_processing = [
        "Variant processing included normalization with VT (Tan et al. 2015)",
        params.run_bcftools ? "and manipulation with BCFtools (Li 2011)." : "."
    ].join(' ').trim()

    def text_annotation = [
        params.run_vep ? "Variant annotation was performed with Ensembl VEP (McLaren et al. 2016)" : "",
        params.run_vcf2maf ? "and conversion to MAF format with vcf2maf (Kandoth et al. 2018)." : "."
    ].join(' ').trim()

    // Combine all sections
    def citation_text = [
        "Tools used in the workflow included:",
        text_preprocessing,
        text_alignment,
        text_alignment_processing,
        text_gatk_preprocessing,
        text_variant_calling,
        text_variant_processing,
        text_annotation,
        "Pipeline results statistics were summarized with MultiQC (Ewels et al. 2016)."
    ].join(' ').trim().replaceAll("\\s+", " ").replaceAll("[,|.] +\\.", ".")

    return citation_text
}


def toolBibliographyText() {
    // Core framework citations (always included)
    def core_citations = [
        "<li>Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. <a href=\"https://doi.org/10.1038/s41587-020-0439-x\">10.1038/s41587-020-0439-x</a></li>",
        "<li>Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. <a href=\"https://doi.org/10.1038/nbt.3820\">10.1038/nbt.3820</a></li>"
    ].join(' ')

    // Quality control and preprocessing
    def qc_preprocessing_citations = [
        "<li>Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. Available: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        params.trimmer == "fastp" ? "<li>Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018 Sep 1;34(17):i884-i890. <a href=\"https://doi.org/10.1093/bioinformatics/bty560\">10.1093/bioinformatics/bty560</a></li>" : ""
    ].join(' ').trim()

    // Alignment tools
    def alignment_citations = [
        params.aligner == "star" ? "<li>Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner Bioinformatics. 2013 Jan 1;29(1):15-21. <a href=\"https://doi.org/10.1093/bioinformatics/bts635\">10.1093/bioinformatics/bts635</a></li>" : "",
        params.aligner == "bwamem" ? "<li>Li H: Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv 2013. <a href=\"https://doi.org/10.48550/arXiv.1303.3997\">10.48550/arXiv.1303.3997</a></li>" : "",
        params.aligner == "bwamem2" ? "<li>M. Vasimuddin, S. Misra, H. Li and S. Aluru, \"Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems,\" 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), 2019, pp. 314-324. <a href=\"https://doi.org/10.1109/IPDPS.2019.00041\">10.1109/IPDPS.2019.00041</a></li>" : "",
        params.aligner == "hisat2" ? "<li>Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019 Aug;37(8):907-915. <a href=\"https://doi.org/10.1038/s41587-019-0201-4\">10.1038/s41587-019-0201-4</a></li>" : ""
    ].join(' ').trim()

    // Alignment processing
    def alignment_processing_citations = [
        "<li>Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. <a href=\"https://doi.org/10.1093/bioinformatics/btp352\">10.1093/bioinformatics/btp352</a></li>",
        params.run_mosdepth ? "<li>Brent S Pedersen, Aaron R Quinlan, Mosdepth: quick coverage calculation for genomes and exomes, Bioinformatics, Volume 34, Issue 5, 01 March 2018, Pages 867–868. <a href=\"https://doi.org/10.1093/bioinformatics/btx699\">10.1093/bioinformatics/btx699</a></li>" : ""
    ].join(' ').trim()

    // GATK preprocessing
    def gatk_citations = [
        params.run_gatk_preprocessing ? "<li>McKenna A, Hanna M, Banks E, et al.: The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. <a href=\"https://doi.org/10.1101/gr.107524.110\">10.1101/gr.107524.110</a></li>" : ""
    ].join(' ').trim()

    // Variant calling
    def variant_calling_citations = [
        params.tools.contains("mutect2") ? "<li>McKenna A, Hanna M, Banks E, et al.: The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010 Sep;20(9):1297-303. <a href=\"https://doi.org/10.1101/gr.107524.110\">10.1101/gr.107524.110</a></li>" : "",
        params.tools.contains("strelka2") ? "<li>Kim S, Scheffler K, Halpern AL, et al.: Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods. 2018 Aug;15(8):591-594. <a href=\"https://doi.org/10.1038/s41592-018-0051-x\">10.1038/s41592-018-0051-x</a></li>" : "",
        params.tools.contains("sage") ? "<li>SAGE is a precise and highly sensitive somatic SNV, MNV and INDEL caller developed by Hartwig Medical Foundation. Available: https://github.com/hartwigmedical/hmftools/tree/master/sage</li>" : ""
    ].join(' ').trim()

    // Variant processing and annotation
    def variant_processing_citations = [
        "<li>Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants. Bioinformatics. 2015 Jul 1;31(13):2202-4. <a href=\"https://doi.org/10.1093/bioinformatics/btv112\">10.1093/bioinformatics/btv112</a></li>",
        params.run_bcftools ? "<li>Li H: A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics. 2011 Nov 1;27(21):2987-93. <a href=\"https://doi.org/10.1093/bioinformatics/btr509\">10.1093/bioinformatics/btr509</a></li>" : "",
        params.run_vep ? "<li>McLaren W, Gil L, Hunt SE, et al.: The Ensembl Variant Effect Predictor. Genome Biol. 2016 Jun 6;17(1):122. <a href=\"https://doi.org/10.1186/s13059-016-0974-4\">10.1186/s13059-016-0974-4</a></li>" : "",
        params.run_vcf2maf ? "<li>Kandoth C, Gao J, Qwangmsk, Mattioni M, Struck A, Boursin Y, Penson A, Chavan S (2018) mskcc/vcf2maf: vcf2maf v1.6.16 (v1.6.16). Zenodo. <a href=\"https://doi.org/10.5281/zenodo.593251\">10.5281/zenodo.593251</a></li>" : ""
    ].join(' ').trim()

    // MultiQC (always included)
    def multiqc_citation = "<li>Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. <a href=\"https://doi.org/10.1093/bioinformatics/btw354\">10.1093/bioinformatics/btw354</a></li>"

    // Software packaging tools (always included)
    def packaging_citations = [
        "<li>Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.</li>",
        "<li>Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. <a href=\"https://doi.org/10.1038/s41592-018-0046-7\">10.1038/s41592-018-0046-7</a></li>",
        "<li>da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. <a href=\"https://doi.org/10.1093/bioinformatics/btx192\">10.1093/bioinformatics/btx192</a></li>",
        "<li>Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2.</li>",
        "<li>Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. <a href=\"https://doi.org/10.1371/journal.pone.0177459\">10.1371/journal.pone.0177459</a></li>"
    ].join(' ')

    // Combine all citations
    def reference_text = [
        core_citations,
        qc_preprocessing_citations,
        alignment_citations,
        alignment_processing_citations,
        gatk_citations,
        variant_calling_citations,
        variant_processing_citations,
        multiqc_citation,
        packaging_citations
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
    meta["tool_citations"]    = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


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
