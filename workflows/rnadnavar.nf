/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    params.multiqc_config,
    params.vep_cache,
    params.star_index,
    params.hisat2_index,
    params.whitelist
    ]

// Validate input parameters
WorkflowRnadnavar.initialise(params, log)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

for (param in checkPathParamList) if (param) file(param, checkIfExists: true)


// Set input, can either be from --input or from automatic retrieval in lib/WorkflowRnadnavar.groovy
if (params.input) {
    ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet("input")
} else {
    ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet("input_restart")
}
// Format samplesheet channel
input_sample = ch_from_samplesheet
        .map{ meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller, maf ->
            // generate patient_sample key to group lanes together
            [ meta.patient + meta.sample, [meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller, maf] ]
        }
        .tap{ ch_with_patient_sample } // save the channel
        .groupTuple() //group by patient_sample to get all lanes
        .map { patient_sample, ch_items ->
            // get number of lanes per sample
            [ patient_sample, ch_items.size() ]
        }
        .combine(ch_with_patient_sample, by: 0) // for each entry add numLanes
        .map { patient_sample, num_lanes, ch_items ->

            (meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller, maf) = ch_items
            if (meta.lane && fastq_2) {
                meta           = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN         = params.seq_center ? "CN:${params.seq_center}\\t" : ''

                def flowcell   = flowcellLaneFromFastq(fastq_1)
                // Don't use a random element for ID, it breaks resuming
                def read_group = "\"@RG\\tID:${flowcell}.${meta.sample}.${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
				if (meta.status >= 2) { // STAR does not need '@RG'
                    read_group  = "ID:${flowcell}.${meta.sample}.${meta.lane} ${CN}PU:${meta.lane} SM:${meta.patient}_${meta.sample} LB:${meta.sample} DS:${params.fasta} PL:${params.seq_platform}"
				}
                meta           = meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'fastq', size: 1]

                if (params.step == 'mapping') return [ meta, [ fastq_1, fastq_2 ] ]
                else {
                    error("Samplesheet contains fastq files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }
            // start for second run
			} else if ((maf || vcf) && params.step=="second_run"){
				if (meta.lane == null) meta.lane = "LX"
				meta            = meta + [id: "${meta.sample}-${meta.lane}-realign".toString()]
                def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}_realign\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
				if (meta.status >= 2) { // STAR does not need '@RG'
					read_group  = "ID:${meta.sample}_${meta.lane}_realign ${CN}PU:${meta.lane} SM:${meta.patient}_${meta.sample} LB:${meta.sample} DS:${params.fasta} PL:${params.seq_platform}"
				}
				if (meta.status >= 2 || meta.status==0){ // these are the files that will go through realignment
	                if (cram)  return [ meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'cram', size: 1], cram, crai, maf ]
	                else if (bam) return [ meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1], bam, bai, maf ]
	                else {
	                    error("Combination error")}
                } else if (meta.status == 1){

                    return [meta + [data_type: 'maf', variantcaller: variantcaller ?: ''], maf]

                }


            // start from BAM
            } else if (meta.lane && bam) {
                if (params.step != 'mapping' && !bai) {
                    error("BAM index (bai) should be provided.")
                }
                meta            = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
				if (meta.status >= 2) { // STAR does not need '@RG'
					read_group  = "ID:${meta.sample}_${meta.lane} ${CN}PU:${meta.lane} SM:${meta.patient}_${meta.sample} LB:${meta.sample} DS:${params.fasta} PL:${params.seq_platform}"
				}
                meta            = meta + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1]

                if (params.step != 'annotate') return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // recalibration
            } else if (table && cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai, table ]
                else {
                    error("Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // recalibration when skipping MarkDuplicates
            } else if (table && bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai, table ]
                else {
                    error("Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // prepare_recalibration or variant_calling
            } else if (cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai ]
                else {
                    error("Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
            } else if (bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']
                if (!(params.step == 'mapping' || params.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is 2 `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }

            // annotation
            } else if (vcf) {
                meta = meta + [id: meta.sample, data_type: 'vcf', variantcaller: variantcaller ?: '']

                if (params.step == 'annotate') return [ meta - meta.subMap('lane'), vcf ]
                else {
                    error("Samplesheet contains vcf files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations")
                }
            } else {
                error("Missing or unknown field in csv file header. Please check your samplesheet")
            }
        }



// Check params logic
if (params.step != 'annotate' && params.tools && !params.build_only_index) {
    // Two checks for ensuring that the pipeline stops with a meaningful error message if
    // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
    // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
    input_sample.filter{ it[0].status == 1 }.ifEmpty{ // In this case, the sample-sheet contains no tumor-samples
        if (!params.build_only_index) {
            def tools_tumor = ['sage','mutect2', 'strelka', 'freebayes']
            def tools_tumor_asked = []
            tools_tumor.each{ tool ->
                if (params.tools.split(',').contains(tool)) tools_tumor_asked.add(tool)
            }
            if (!tools_tumor_asked.isEmpty()) {
                error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
            }
        }
    }
    input_sample.filter{ it[0].status == 0 }.ifEmpty{ // In this case, the sample-sheet contains no normal/germline-samples
        def tools_requiring_normal_samples = ['sage','mutect2', 'strelka', 'freebayes'] // Will implement tumour only in the near future
        def requested_tools_requiring_normal_samples = []
        tools_requiring_normal_samples.each{ tool_requiring_normal_samples ->
            if (params.tools.split(',').contains(tool_requiring_normal_samples)) requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
        }
        if (!requested_tools_requiring_normal_samples.isEmpty()) {
            error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample : ' + requested_tools_requiring_normal_samples.join(", "))
        }
    }
}

// Fails when wrongful extension for intervals file
if (params.wes && !params.step == 'annotate') {
    if (params.intervals && !params.intervals.endsWith("bed"))  error("Target file specified with `--intervals` must be in BED format for targeted data")
    else log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
} else if (params.intervals && !params.intervals.endsWith("bed") && !params.intervals.endsWith("list")) error("Intervals file must end with .bed, .list, or .interval_list")


// Fails when missing params for STAR
if (!params.star_index && !params.gtf && !params.gff){
     exit 1,"GTF|GFF3 file is required to build a STAR reference index! Use option --gtf|--gff to provide a GTF|GFF file."
}

// Warns when missing files or params for mutect2
if(params.tools && params.tools.split(',').contains('mutect2')){
    if(!params.pon){
        log.warn("No Panel-of-normal was specified for Mutect2.\nIt is highly recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2\nFor more information on how to create one: https://gatk.broadinstitute.org/hc/en-us/articles/5358921041947-CreateSomaticPanelOfNormals-BETA-")
    }
    if(!params.germline_resource){
        log.warn("If Mutect2 is specified without a germline resource, no filtering will be done.\nIt is recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2")
    }
    if(params.pon && params.pon.contains("/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz")){
        log.warn("The default Panel-of-Normals provided by GATK is used for Mutect2.\nIt is highly recommended to generate one from normal samples that are technical similar to the tumor ones.\nFor more information: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-")
    }
}

if(params.tools && params.tools.split(',').contains('sage')){
	if(!params.sage_ensembl_dir){
        log.error "SAGE requires ensembl resource file. Please provide `--sage_ensembl_dir`\nYou can skip this step in the workflow by removing sage from `--tools` to the command."
        exit 1
	}
}


// Fails when missing resources for baserecalibrator
if(!params.dbsnp && !params.known_indels){
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!params.skip_tools || (params.skip_tools && !params.skip_tools.split(',').contains('baserecalibrator')))){
        log.error "Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command."
        exit 1
    }
}

// Fails when missing tools for variant_calling or annotate
if ((params.step == 'variant_calling' || params.step == 'annotate') && !params.tools) {
    error("Please specify at least one tool when using `--step ${params.step}`.\nhttps://nf-co.re/rnadnavar/parameters#tools")
}

if ((params.download_cache) && (params.snpeff_cache || params.vep_cache)) {
    error("Please specify either `--download_cache` or `--vep_cache`.\nhttps://nf-co.re/rnadnavar/dev/usage#how-to-customise-vep-annotation")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// Build the genome index and other reference files
include { PREPARE_REFERENCE_AND_INTERVALS   } from '../subworkflows/local/prepare_reference_and_intervals/main'
// Download annotation cache if needed
include { ENSEMBLVEP_DOWNLOAD               } from '../modules/nf-core/ensemblvep/download/main'

// Alignment
include { BAM_ALIGN                         } from '../subworkflows/local/bam_align/main'

// Core subworkflows of the pipeline
include { BAM_VARIANT_CALLING_PRE_POST_PROCESSING               } from '../subworkflows/local/bam_variant_calling_pre_post_processing/main'

// Second run
include { BAM_EXTRACT_READS_HISAT2_ALIGN as PREPARE_SECOND_RUN  } from '../subworkflows/local/prepare_second_run/main'
include { BAM_VARIANT_CALLING_PRE_POST_PROCESSING as SECOND_RUN } from '../subworkflows/local/bam_variant_calling_pre_post_processing/main'
//include { FILTERING_RNA          } from '../subworkflows/local/rna_filtering'
//
//
// MODULE: Installed directly from nf-core/modules
//
//FASTQC
include { FASTQC                                      } from '../modules/nf-core/fastqc/main'
// MULTIQC
include { MULTIQC                                     } from '../modules/nf-core/multiqc/main'
// REPORTING VERSIONS OF SOFTWARE USED
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/custom/dumpsoftwareversions/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
========================================================================================
   VARIABLES
========================================================================================
*/

whitelist  = params.whitelist  ? Channel.fromPath(params.whitelist).collect() : Channel.value([])
blacklist  = params.blacklist  ? Channel.fromPath(params.blacklist).collect() : Channel.value([])
darned     = params.darned     ? Channel.fromPath(params.darned).collect()    : Channel.value([])
radar      = params.radar      ? Channel.fromPath(params.radar).collect()     : Channel.value([])
nat        = params.nat        ? Channel.fromPath(params.nat).collect()       : Channel.value([])
redi       = params.redi       ? Channel.fromPath(params.redi).collect()      : Channel.value([])

// Info required for completion email and summary
def multiqc_report = []

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNADNAVAR {

    // Initialise MULTIQC
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    // To gather all QC reports for MultiQC
    reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

    // Download cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache    ? [] : Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])

    if (params.download_cache) {
        ENSEMBLVEP_DOWNLOAD(ensemblvep_info)
        vep_cache    = ENSEMBLVEP_DOWNLOAD.out.cache.collect().map{ meta, cache -> [ cache ] }

        versions = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
    }

	// STEP 0: Build reference and indices if needed
	PREPARE_REFERENCE_AND_INTERVALS()
	versions = versions.mix(PREPARE_REFERENCE_AND_INTERVALS.out.versions)

	// Reference and intervals variables
	fasta                         =    PREPARE_REFERENCE_AND_INTERVALS.out.fasta
	fasta_fai                     =    PREPARE_REFERENCE_AND_INTERVALS.out.fasta_fai
	dict                          =    PREPARE_REFERENCE_AND_INTERVALS.out.dict
	germline_resource             =    PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource
	germline_resource_tbi         =    PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource_tbi
	intervals                     =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals
	intervals_for_preprocessing   =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals_for_preprocessing
	// specific for variant calling
	intervals_bed_combined        =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_combined
	intervals_bed_gz_tbi          =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_gz_tbi
	intervals_bed_gz_tbi_combined =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_gz_tbi_combined
	dbsnp                         =    PREPARE_REFERENCE_AND_INTERVALS.out.dbsnp
	dbsnp_tbi                     =    PREPARE_REFERENCE_AND_INTERVALS.out.dbsnp_tbi
	pon                           =    PREPARE_REFERENCE_AND_INTERVALS.out.pon
	pon_tbi                       =    PREPARE_REFERENCE_AND_INTERVALS.out.pon_tbi
	known_sites_indels            =    PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_indels
	known_sites_indels_tbi        =    PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_indels_tbi
	known_sites_snps              =    PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_snps
	known_sites_snps_tbi          =    PREPARE_REFERENCE_AND_INTERVALS.out.known_sites_snps_tbi

	intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }
// STEP 1: ALIGNMENT PREPROCESSING
	BAM_ALIGN(
	   PREPARE_REFERENCE_AND_INTERVALS.out.bwa,
	   PREPARE_REFERENCE_AND_INTERVALS.out.bwamem2,
	   PREPARE_REFERENCE_AND_INTERVALS.out.dragmap,
	   PREPARE_REFERENCE_AND_INTERVALS.out.star_index,
	   PREPARE_REFERENCE_AND_INTERVALS.out.gtf,
	   input_sample
	   )

	reports = reports.mix(BAM_ALIGN.out.reports)
	versions = versions.mix(BAM_ALIGN.out.versions)
	// 5 MAIN STEPS: GATK PREPROCESING - VARIANT CALLING - NORMALIZATION - CONSENSUS - ANNOTATION
	BAM_VARIANT_CALLING_PRE_POST_PROCESSING(
	    input_sample,              // input from CSV if applicable
	    BAM_ALIGN.out.bam_mapped,  // input from mapping
	    BAM_ALIGN.out.cram_mapped,  // input from mapping
	    fasta,                     // fasta reference file
	    fasta_fai,                 // fai for fasta file
	    dict,                      //
	    dbsnp,
	    dbsnp_tbi,
	    pon,
	    pon_tbi,
	    known_sites_indels,
	    known_sites_indels_tbi,
	    germline_resource,
	    germline_resource_tbi,
	    intervals,
	    intervals_for_preprocessing,
	    intervals_bed_gz_tbi,
	    intervals_bed_combined,
	    intervals_and_num_intervals,
	    intervals_bed_gz_tbi_combined,
	    null,  // to repeat rescue consensus TODO: is this the best strategy?
	    null,  // to repeat rescue consensus
	    false  // is second run
	)

    reports  = reports.mix(BAM_VARIANT_CALLING_PRE_POST_PROCESSING.out.reports)
    versions = versions.mix(BAM_VARIANT_CALLING_PRE_POST_PROCESSING.out.versions)


    version_yaml = Channel.empty()
    if (!(params.skip_tools && params.skip_tools.split(',').contains('versions'))) {
        CUSTOM_DUMPSOFTWAREVERSIONS(versions.unique().collectFile(name: 'collated_versions.yml'))
        version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }

    if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {
        workflow_summary    = WorkflowRnadnavar.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        methods_description    = WorkflowRnadnavar.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
        ch_methods_description = Channel.value(methods_description)

        multiqc_files = Channel.empty()
        multiqc_files = multiqc_files.mix(version_yaml)
        multiqc_files = multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        multiqc_files = multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        multiqc_files = multiqc_files.mix(reports.collect().ifEmpty([]))

        MULTIQC(multiqc_files.collect(), ch_multiqc_config.collect().ifEmpty([]), ch_multiqc_custom_config.collect().ifEmpty([]), ch_multiqc_logo.collect().ifEmpty([]))

        multiqc_report = MULTIQC.out.report.toList()
        versions = versions.mix(MULTIQC.out.versions)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
