/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnadnavar.initialise(params, log)

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
    params.snpeff_cache,
    params.vep_cache,
    params.star_index,
    params.hisat2_index,
    params.whitelist
    ]

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
for (param in checkPathParamList) {
    if (param) {
        file(param, checkIfExists: true)
        }
    }

// Set input, can either be from --input or from automatic retrieval in lib/WorkflowRnadnavar.groovy
ch_input_sample = extract_csv(file(params.input))

// Fails when wrongful extension for intervals file
if (params.wes && !params.step == 'annotate') {
    if (params.intervals && !params.intervals.endsWith("bed")) exit 1, "Target file specified with `--intervals` must be in BED format for targeted data"
    else log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
} else if (params.intervals && !params.intervals.endsWith("bed") && !params.intervals.endsWith("interval_list")) exit 1, "Intervals file must end with .bed or .interval_list"

if(params.step == 'mapping' && params.aligner.contains("dragmap") && !(params.skip_tools && params.skip_tools.split(',').contains("baserecalibrator"))){
    log.warn("DragMap was specified as aligner. Base recalibration is not contained in --skip_tools. It is recommended to skip baserecalibration when using DragMap\nhttps://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode")
}

// Fails when missing params for STAR
if (!params.star_index && !params.gtf && !params.gff)
    {
    exit 1,
    "GTF|GFF3 file is required to build a STAR reference index! Use option --gtf|--gff to provide a GTF|GFF file."
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

// Fails when missing resources for baserecalibrator
if(!params.dbsnp && !params.known_indels){
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!params.skip_tools || (params.skip_tools && !params.skip_tools.split(',').contains('baserecalibrator')))){
        log.error "Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command."
        exit 1
    }
}

// Fails when missing tools for variant_calling or annotate
if ((params.step == 'variant_calling' || params.step == 'annotate') && !params.tools) {
    log.error "Please specify at least one tool when using `--step ${params.step}`.\nhttps://nf-co.re/rnadnavar/parameters#tools"
    exit 1
}

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[params.genome]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

file("${params.outdir}").mkdirs()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = [
                            file("$projectDir/assets/multiqc_config.yml", checkIfExists: true),
                            file("$projectDir/assets/nf-core-rnadnavar_logo_light.png", checkIfExists: true)
                            ]
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CORE_RUN                                            } from '../subworkflows/local/core_workflow_pass'
include { CORE_RUN as SECOND_RUN                              } from '../subworkflows/local/core_workflow_pass'

// Build the genome index and other reference files
include { PREPARE_REFERENCE_AND_INTERVALS                     } from '../subworkflows/local/prepare_reference_and_intervals'
include { MAPPING                                             } from '../subworkflows/local/mapping'

// Filtering
include { PREPARE_SECOND_RUN     } from '../subworkflows/local/prepare_second_run'
include { FILTERING_RNA          } from '../subworkflows/local/rna_filtering'
// REPORTING VERSIONS OF SOFTWARE USED
include { CUSTOM_DUMPSOFTWAREVERSIONS                          } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

// MULTIQC
include { MULTIQC                                              } from '../modules/nf-core/modules/multiqc/main'




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

    // To gather all QC reports for MultiQC
    ch_reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    ch_versions = Channel.empty()

// STEP 0: Build reference and indices if needed
    PREPARE_REFERENCE_AND_INTERVALS()
    ch_versions = ch_versions.mix(PREPARE_REFERENCE_AND_INTERVALS.out.versions)

    // Reference and intervals variables
    fasta                       =    PREPARE_REFERENCE_AND_INTERVALS.out.fasta
    fasta_fai                   =    PREPARE_REFERENCE_AND_INTERVALS.out.fasta_fai
    dict                        =    PREPARE_REFERENCE_AND_INTERVALS.out.dict
    germline_resource           =    PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource
    germline_resource_tbi       =    PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource_tbi
    intervals                   =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals
    intervals_for_preprocessing =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals_for_preprocessing
    ch_interval_list_split      =    PREPARE_REFERENCE_AND_INTERVALS.out.ch_interval_list_split
    // specific for variant calling
    intervals_bed_combined      =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi        =    PREPARE_REFERENCE_AND_INTERVALS.out.intervals_bed_gz_tbi
    dbsnp                       =    PREPARE_REFERENCE_AND_INTERVALS.out.dbsnp
    dbsnp_tbi                   =    PREPARE_REFERENCE_AND_INTERVALS.out.dbsnp_tbi
    pon                         =    PREPARE_REFERENCE_AND_INTERVALS.out.pon
    pon_tbi                     =    PREPARE_REFERENCE_AND_INTERVALS.out.pon_tbi
    germline_resource           =    PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource
    germline_resource_tbi       =    PREPARE_REFERENCE_AND_INTERVALS.out.germline_resource_tbi


// STEP 1: ALIGNMENT PREPROCESSING
    MAPPING(
        PREPARE_REFERENCE_AND_INTERVALS.out.bwa,
        PREPARE_REFERENCE_AND_INTERVALS.out.bwamem2,
        PREPARE_REFERENCE_AND_INTERVALS.out.dragmap,
        PREPARE_REFERENCE_AND_INTERVALS.out.star_index,
        PREPARE_REFERENCE_AND_INTERVALS.out.gtf,
        ch_input_sample
        )
    ch_reports = ch_reports.mix(MAPPING.out.reports)
    ch_versions = ch_versions.mix(MAPPING.out.versions)

    // 5 MAIN STEPS: GATK PREPROCESING - VARIANT CALLING - NORMALIZATION - CONSENSUS - ANNOTATION
    CORE_RUN(
        params.step,
        params.skip_tools,
        ch_input_sample,           // input from CSV if applicable
        MAPPING.out.ch_bam_mapped, // input from mapping
        fasta,                     // fasta reference file
        fasta_fai,                 // fai for fasta file
        dict,                      //
        dbsnp,
        dbsnp_tbi,
        pon,
        pon_tbi,
        germline_resource,
        germline_resource_tbi,
        intervals,
        intervals_for_preprocessing,
        ch_interval_list_split,
        intervals_bed_gz_tbi,
        intervals_bed_combined,
        null,  // to repeat rescue consensus
        null  // to repeat rescue consensus
    )

    ch_reports = ch_reports.mix(CORE_RUN.out.reports)
    ch_versions = ch_versions.mix(CORE_RUN.out.versions)

    PREPARE_SECOND_RUN(ch_input_sample,           // input from CSV if applicable
                           CORE_RUN.out.maf,
                           MAPPING.out.bwa_bams, // for dna re-alignments
                           MAPPING.out.star_bams,  // for rnare-alignments
                           fasta,
                           fasta_fai,
                           dict,
                           PREPARE_REFERENCE_AND_INTERVALS.out.hisat2_index,
                           PREPARE_REFERENCE_AND_INTERVALS.out.splicesites
                           ) // do mapping with hisat2

    ch_reports = ch_reports.mix(PREPARE_SECOND_RUN.out.reports)
    ch_versions = ch_versions.mix(PREPARE_SECOND_RUN.out.versions)

    SECOND_RUN(
        "markduplicates",                      // step to start with
        "baserecalibrator,baserecalibrator_report,contamination,learnreadorientation",
        ch_input_sample,                       // input from CSV if applicable
        PREPARE_SECOND_RUN.out.ch_bam_mapped, // input from mapping
        fasta,                                 // fasta reference file
        fasta_fai,                             // fai for fasta file
        dict,                                    //
        dbsnp,
        dbsnp_tbi,
        pon,
        pon_tbi,
        germline_resource,
        germline_resource_tbi,
        intervals,
        intervals_for_preprocessing,
        ch_interval_list_split,
        intervals_bed_gz_tbi,
        intervals_bed_combined,
        CORE_RUN.out.vcf_consensus_dna,  // to repeat rescue consensus
        CORE_RUN.out.vcfs_status_dna  // to repeat rescue consensus
    )

    ch_reports = ch_reports.mix(SECOND_RUN.out.reports)
    ch_versions = ch_versions.mix(SECOND_RUN.out.versions)

    FILTERING_RNA(CORE_RUN.out.maf_rna,
                  SECOND_RUN.out.maf_rna,
                  fasta)
    ch_versions = ch_versions.mix(FILTERING_RNA.out.versions)

// REPORTING

    ch_version_yaml = Channel.empty()
    if (!(params.skip_tools && params.skip_tools.split(',').contains('versions'))) {
        CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
        ch_version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }

    // MODULE: MultiQC
    // Present summary of reads, alignment, duplicates, BSQR stats for all samples as well as workflow summary/parameters as single report
    if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {
        workflow_summary    = WorkflowRnadnavar.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files =  Channel.empty().mix(ch_version_yaml,
                                            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                                            ch_reports.collect().ifEmpty([]))

        ch_multiqc_configs = Channel.from(ch_multiqc_config).mix(ch_multiqc_custom_config).ifEmpty([])

        MULTIQC(ch_multiqc_files.collect(), ch_multiqc_configs.collect())
        multiqc_report = MULTIQC.out.report.toList()
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Function to extract information (meta data + file(s)) from csv file(s)
def extract_csv(csv_file) {

    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {numberOfLinesInSampleSheet++}
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Samplesheet had less than two lines. The sample sheet must be a csv file with a header, so at least two lines."
            System.exit(1)
        }
    }

    // Additional check of sample sheet:
    // 1. If params.step == "mapping", then each row should specify a lane and the same combination of patient, sample and lane shouldn't be present in different rows.
    // 2. The same sample shouldn't be listed for different patients.
    def patient_sample_lane_combinations_in_samplesheet = []
    def sample2patient = [:]

    Channel.from(csv_file).splitCsv(header: true)
        .map{ row ->
            if (params.step == "mapping") {
                if ( !row.lane ) {  // This also handles the case where the lane is left as an empty string
                    log.error('The sample sheet should specify a lane for patient "' + row.patient.toString() + '" and sample "' + row.sample.toString() + '".')
                    System.exit(1)
                }
                def patient_sample_lane = [row.patient.toString(), row.sample.toString(), row.lane.toString()]
                if (patient_sample_lane in patient_sample_lane_combinations_in_samplesheet) {
                    log.error('The patient-sample-lane combination "' + row.patient.toString() + '", "' + row.sample.toString() + '", and "' + row.lane.toString() + '" is present multiple times in the sample sheet.')
                    System.exit(1)
                } else {
                    patient_sample_lane_combinations_in_samplesheet.add(patient_sample_lane)
                }
            }
            if (!sample2patient.containsKey(row.sample.toString())) {
                sample2patient[row.sample.toString()] = row.patient.toString()
            } else if (sample2patient[row.sample.toString()] != row.patient.toString()) {
                log.error('The sample "' + row.sample.toString() + '" is registered for both patient "' + row.patient.toString() + '" and "' + sample2patient[row.sample.toString()] + '" in the sample sheet.')
                System.exit(1)
            }
        }
    // keep count of the number of samples
    sample_count_all = 0
    sample_count_normal = 0
    sample_count_tumor = 0
    sample_count_rna = 0

    Channel.from(csv_file).splitCsv(header: true)
        // Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            sample_count_all++
            if (!(row.patient && row.sample)){
                log.error "Missing field in csv file header. The csv file must have fields named 'patient' and 'sample'."
                System.exit(1)
            }
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing

        def meta = [:]

        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.patient) meta.patient = row.patient.toString()
        if (row.sample)  meta.sample  = row.sample.toString()

        // If no status specified, sample is assumed normal
        if (row.status) meta.status = row.status.toInteger()
        else meta.status = 0

        if (meta.status == 0) sample_count_normal++
        else if (meta.status == 1) sample_count_tumor++  // TODO check if elif is valid in here
        else sample_count_rna++
        // TODO: think about what other condition we will have here now
        // Two checks for ensuring that the pipeline stops with a meaningful error message if
        // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
        // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
        if ((sample_count_normal == sample_count_all) && params.tools) { // In this case, the sample-sheet contains no tumor-samples
            def tools_tumor = ['sage', 'mutect2', 'strelka2']  // This will be applied to tumour DNA and tumour RNA
            def tools_tumor_asked = []
            tools_tumor.each{ tool ->
                if (params.tools.split(',').contains(tool)) tools_tumor_asked.add(tool)
            }
            if (!tools_tumor_asked.isEmpty()) {
                log.error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
                System.exit(1)
            }
        // TODO no need to do anything with the germline - can this be removed?
        } else if ((sample_count_tumor == sample_count_all) && params.tools) {  // In this case, the sample-sheet contains no normal/germline-samples
            def tools_requiring_normal_samples = ['ascat', 'deepvariant', 'haplotypecaller']
            def requested_tools_requiring_normal_samples = []
            tools_requiring_normal_samples.each{ tool_requiring_normal_samples ->
                if (params.tools.split(',').contains(tool_requiring_normal_samples)) requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
            }
            if (!requested_tools_requiring_normal_samples.isEmpty()) {
                log.error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample : ' + requested_tools_requiring_normal_samples.join(", "))
                System.exit(1)
            }
        }

        // mapping with fastq
        if (row.lane && row.fastq_2) {
            meta.id         = "${row.sample}-${row.lane}".toString()
            def fastq_1     = file(row.fastq_1, checkIfExists: true)
            def fastq_2     = file(row.fastq_2, checkIfExists: true)
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''

            def flowcell    = flowcellLaneFromFastq(fastq_1)
            //Don't use a random element for ID, it breaks resuming
            def read_group  = "\"@RG\\tID:${flowcell}.${row.sample}.${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.patient}_${row.sample}\\tLB:${row.sample}\\tDS:${params.fasta}\\tPL:${params.seq_platform}\""
            if (meta.status == 2) { // STAR does not need '@RG'
               read_group  = "ID:${flowcell}.${row.sample}.${row.lane} ${CN}PU:${row.lane} SM:${row.patient}_${row.sample} LB:${row.sample} DS:${params.fasta} PL:${params.seq_platform}"
            }

            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'fastq'

            meta.size       = 1 // default number of splitted fastq

            if (params.step == 'mapping') return [meta, [fastq_1, fastq_2]]
            else {
                log.error "Samplesheet contains fastq files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // start from BAM
        } else if (row.lane && row.bam) {
            if (!row.bai) {
                log.error "BAM index (bai) should be provided."
            }
            meta.id         = "${row.sample}-${row.lane}".toString()
            def bam         = file(row.bam,   checkIfExists: true)
            def bai         = file(row.bai,   checkIfExists: true)
            def CN          = params.seq_center ? "CN:${params.seq_center}\\t" : ''
            def read_group  = "\"@RG\\tID:${row.sample}_${row.lane}\\t${CN}PU:${row.lane}\\tSM:${row.sample}\\tLB:${row.sample}\\tPL:${params.seq_platform}\""
            if (meta.status == 2) { // STAR does not need '@RG'
               read_group  = "ID:${row.sample}_${row.lane} ${CN}PU:${row.lane} SM:${row.sample} LB:${row.sample} PL:${params.seq_platform}"
            }

            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'bam'

            meta.size       = 1 // default number of splitted fastq

            if (params.step != 'annotate') return [meta, bam, bai]
            else {
                log.error "Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // recalibration
        } else if (row.table && row.cram) {
            meta.id   = meta.sample
            def cram  = file(row.cram,  checkIfExists: true)
            def crai  = file(row.crai,  checkIfExists: true)
            def table = file(row.table, checkIfExists: true)

            meta.data_type  = 'cram'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, cram, crai, table]
            else {
                log.error "Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // recalibration when skipping MarkDuplicates
        } else if (row.table && row.bam) {
            meta.id   = meta.sample
            def bam   = file(row.bam,   checkIfExists: true)
            def bai   = file(row.bai,   checkIfExists: true)
            def table = file(row.table, checkIfExists: true)

            meta.data_type  = 'bam'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, bam, bai, table]
            else {
                log.error "Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // prepare_recalibration or variant_calling
        } else if (row.cram) {
            meta.id = meta.sample
            def cram = file(row.cram, checkIfExists: true)
            def crai = file(row.crai, checkIfExists: true)

            meta.data_type  = 'cram'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, cram, crai]
            else {
                log.error "Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
        } else if (row.bam) {
            meta.id = meta.sample
            def bam = file(row.bam, checkIfExists: true)
            def bai = file(row.bai, checkIfExists: true)

            meta.data_type  = 'bam'

            if (!(params.step == 'mapping' || params.step == 'annotate')) return [meta, bam, bai]
            else {
                log.error "Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations"
                System.exit(1)
            }

        // annotation
        } else if (row.vcf) {
            meta.id = meta.sample
            def vcf = file(row.vcf, checkIfExists: true)

            meta.data_type     = 'vcf'
            meta.variantcaller = row.variantcaller ?: ''

            if (params.step == 'annotate') return [meta, vcf]
            else {
                log.error "Samplesheet contains vcf files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/rnadnavar/usage#input-samplesheet-configurations"
                System.exit(1)
            }
        } else {
            log.error "Missing or unknown field in csv file header. Please check your samplesheet"
            System.exit(1)
        }
    }
}
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
