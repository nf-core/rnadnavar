/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnadnavar.initialise(params, log)

// TODO nf-core: For now added paramers from rnavar to start with
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.fasta,
    params.fasta_fai,
    params.dict,
    params.bwa,
    params.bwamem2,
    params.gtf,
    params.gff,
    params.dbsnp,
    params.dbsnp_tbi,
    params.known_indels,
    params.known_indels_tbi,
    params.snpeff_cache,
    params.vep_cache,
    params.star_index
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
ch_input_sample = extract_csv(file(params.input, checkIfExists: true))

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
// Warns when missing resources for haplotypecaller
if(!params.dbsnp && !params.known_indels){
    if (params.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!params.skip_tools || (params.skip_tools && !params.skip_tools.split(',').contains('baserecalibrator')))){
        log.error "Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command."
        exit 1
    }
    if(params.tools && params.tools.split(',').contains('haplotypecaller')){
        log.warn "If Haplotypecaller is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-"

    }
}
if (params.joint_germline && (!params.dbsnp || !params.known_indels || !params.known_snps || params.no_intervals)){
    log.warn "If Haplotypecaller is specified, without `--dbsnp`, `--known_snps`, `--known_indels` or the associated resource labels (ie `known_snps_vqsr`), no variant recalibration will be done. For recalibration you must provide all of these resources.\nFor more information see VariantRecalibration: https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator \nJoint germline variant calling also requires intervals in order to genotype the samples. As a result, if `--no_intervals` is set to `true` the joint germline variant calling will not be performed."
}

// Fails when missing tools for variant_calling or annotate
if ((params.step == 'variant_calling' || params.step == 'annotate') && !params.tools) {
    log.error "Please specify at least one tool when using `--step ${params.step}`.\nhttps://nf-co.re/sarek/parameters#tools"
    exit 1
}

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[params.genome]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
ch_rnadnavar_logo           = Channel.fromPath(file("$projectDir/assets/nf-core-rnadnavar_logo_light.png", checkIfExists: true))


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
dbsnp              = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()                    : Channel.value([])
known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()               : Channel.value([])
fasta              = params.fasta              ? Channel.fromPath(params.fasta).collect()                    : Channel.empty()
fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()                : Channel.empty()
germline_resource  = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect()        : Channel.value([]) //Mutec2 does not require a germline resource, so set to optional input
known_indels       = params.known_indels       ? Channel.fromPath(params.known_indels).collect()             : Channel.value([])
known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()               : Channel.value([])
pon                = params.pon                ? Channel.fromPath(params.pon).collect()                      : Channel.value([]) //PON is optional for Mutect2 (but highly recommended)


// Create samplesheets to restart from different steps
include { MAPPING_CSV                                          } from '../subworkflows/local/mapping_csv'
include { MARKDUPLICATES_CSV                                   } from '../subworkflows/local/markduplicates_csv'
include { PREPARE_RECALIBRATION_CSV                            } from '../subworkflows/local/prepare_recalibration_csv'
include { RECALIBRATE_CSV                                      } from '../subworkflows/local/recalibrate_csv'
include { VARIANTCALLING_CSV                                   } from '../subworkflows/local/variantcalling_csv'

// Build the genome index and other reference files
include { PREPARE_GENOME                                       } from '../subworkflows/local/prepare_genome'

// Build intervals if needed
include { PREPARE_INTERVALS                                    } from '../subworkflows/local/prepare_intervals'

// Convert BAM files to FASTQ files
include { ALIGNMENT_TO_FASTQ as ALIGNMENT_TO_FASTQ_INPUT       } from '../subworkflows/nf-core/alignment_to_fastq'
include { ALIGNMENT_TO_FASTQ as ALIGNMENT_TO_FASTQ_UMI         } from '../subworkflows/nf-core/alignment_to_fastq'

// Run FASTQC
include { RUN_FASTQC                                           } from '../subworkflows/nf-core/run_fastqc'

// TRIM/SPLIT FASTQ Files
include { FASTP                                                } from '../modules/nf-core/modules/fastp/main'

// Create umi consensus bams from fastq
include { CREATE_UMI_CONSENSUS                                 } from '../subworkflows/nf-core/fgbio_create_umi_consensus/main'

// Map input reads to reference genome
include { GATK4_MAPPING                                        } from '../subworkflows/nf-core/gatk4/mapping/main'

// Merge and index BAM files (optional)
include { MERGE_INDEX_BAM                                      } from '../subworkflows/nf-core/merge_index_bam'

include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM               } from '../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM_RECAL         } from '../modules/nf-core/modules/samtools/convert/main'

include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM               } from '../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM_VARIANTCALLING} from '../modules/nf-core/modules/samtools/convert/main'

// Mark Duplicates (+QC)
include { MARKDUPLICATES                                       } from '../subworkflows/nf-core/gatk4/markduplicates/main'

// Mark Duplicates SPARK (+QC)
include { MARKDUPLICATES_SPARK                                 } from '../subworkflows/nf-core/gatk4/markduplicates_spark/main'

// Convert to CRAM (+QC)
include { BAM_TO_CRAM                                          } from '../subworkflows/nf-core/bam_to_cram'

// QC on CRAM
include { CRAM_QC                                              } from '../subworkflows/nf-core/cram_qc'

// Create recalibration tables
include { PREPARE_RECALIBRATION                                } from '../subworkflows/nf-core/gatk4/prepare_recalibration/main'

// Create recalibration tables SPARK
include { PREPARE_RECALIBRATION_SPARK                          } from '../subworkflows/nf-core/gatk4/prepare_recalibration_spark/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { RECALIBRATE                                          } from '../subworkflows/nf-core/gatk4/recalibrate/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { RECALIBRATE_SPARK                                    } from '../subworkflows/nf-core/gatk4/recalibrate_spark/main'

// Variant calling on tumor/normal pair TODO: add tumour only for next version
include { PAIR_VARIANT_CALLING                                 } from '../subworkflows/local/pair_variant_calling'

include { VCF_QC                                               } from '../subworkflows/nf-core/vcf_qc'

// Annotation
include { ANNOTATE                                             } from '../subworkflows/local/annotate'

// REPORTING VERSIONS OF SOFTWARE USED
include { CUSTOM_DUMPSOFTWAREVERSIONS                          } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

// MULTIQC
include { MULTIQC                                              } from '../modules/nf-core/modules/multiqc/main'


/*
========================================================================================
    IMPORT NF-CORE SUBWORKFLOWS
========================================================================================
*/

include { ALIGN_STAR                    } from '../subworkflows/nf-core/align_star'         // Align reads to genome and sort and index the alignment file
include { SPLITNCIGAR                   } from '../subworkflows/nf-core/splitncigar'        // Splits reads that contain Ns in their cigar string
include { SPLITNCIGAR2                   } from '../subworkflows/nf-core/splitncigar2'        // Splits reads that contain Ns in their cigar string - without intervals option
// TODO not needed?
// include { RECALIBRATE                   } from '../subworkflows/nf-core/recalibrate'        // Estimate and correct systematic bias

/*
========================================================================================
    VARIABLES
========================================================================================
*/

// Check STAR alignment parameters
def prepareToolIndices  = params.aligner
def seq_platform        = params.seq_platform ? params.seq_platform : []
def seq_center          = params.seq_center ? params.seq_center : []

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

    // Build indices if needed
    PREPARE_GENOME(
        dbsnp,
        fasta,
        fasta_fai,
        germline_resource,
        known_indels,
        known_snps,
        pon)

    // Gather built indices or get them from the params
    bwa                    = params.fasta                   ? params.bwa                        ? Channel.fromPath(params.bwa).collect()                   : PREPARE_GENOME.out.bwa                   : []
    bwamem2                = params.fasta                   ? params.bwamem2                    ? Channel.fromPath(params.bwamem2).collect()               : PREPARE_GENOME.out.bwamem2               : []
    dragmap                = params.fasta                   ? params.dragmap                    ? Channel.fromPath(params.dragmap).collect()               : PREPARE_GENOME.out.hashtable             : []
    dict                   = params.fasta                   ? params.dict                       ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
    fasta_fai              = params.fasta                   ? params.fasta_fai                  ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []
    dbsnp_tbi              = params.dbsnp                   ? params.dbsnp_tbi                  ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.value([])
    germline_resource_tbi  = params.germline_resource       ? params.germline_resource_tbi      ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
    known_indels_tbi       = params.known_indels            ? params.known_indels_tbi           ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.value([])
    known_snps_tbi         = params.known_snps              ? params.known_snps_tbi             ? Channel.fromPath(params.known_snps_tbi).collect()        : PREPARE_GENOME.out.known_snps_tbi        : Channel.value([])
    pon_tbi                = params.pon                     ? params.pon_tbi                    ? Channel.fromPath(params.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : []

    // Gather index for mapping given the chosen aligner
    ch_map_index = params.aligner == "bwa-mem" ? bwa :
        params.aligner == "bwa-mem2" ? bwamem2 :
        dragmap

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    known_sites_snps     = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals_bed_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined  // [interval.bed] all intervals in one file
    intervals_for_preprocessing = params.wes          ? intervals_bed_combined : []      // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)

    intervals                   = PREPARE_INTERVALS.out.intervals_bed        // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi        = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [interval_bed, tbi, num_intervals] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    // Gather used softwares versions
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

    // PREPROCESSING

    if (params.step == 'mapping') {

        // Figure out if input is bam or fastq
        ch_input_sample.branch{
            bam:   it[0].data_type == "bam"
            fastq: it[0].data_type == "fastq"
        }.set{ch_input_sample_type}

        // convert any bam input to fastq
        // Fasta are not needed when converting bam to fastq -> []
        ALIGNMENT_TO_FASTQ_INPUT(ch_input_sample_type.bam, [])

        // gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        ch_input_fastq = ch_input_sample_type.fastq.mix(ALIGNMENT_TO_FASTQ_INPUT.out.reads)

        // STEP 0: QC & TRIM
        // `--skip_tools fastqc` to skip fastqc
        // trim only with `--trim_fastq`
        // additional options to be set up

        // QC
        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            RUN_FASTQC(ch_input_fastq)

            ch_reports  = ch_reports.mix(RUN_FASTQC.out.fastqc_zip.collect{meta, logs -> logs})
            ch_versions = ch_versions.mix(RUN_FASTQC.out.versions)
        }

        // UMI consensus calling
        if (params.umi_read_structure) {
            CREATE_UMI_CONSENSUS(
                ch_input_fastq,
                fasta,
                ch_map_index,
                umi_read_structure,
                params.group_by_umi_strategy
            )

            bamtofastq = CREATE_UMI_CONSENSUS.out.consensusbam.map{meta, bam -> [meta,bam,[]]}

            // convert back to fastq for further preprocessing
            ALIGNMENT_TO_FASTQ_UMI(bamtofastq, [])

            ch_reads_fastp = ALIGNMENT_TO_FASTQ_UMI.out.reads

            // Gather used softwares versions
            ch_versions = ch_versions.mix(ALIGNMENT_TO_FASTQ_UMI.out.versions)
            ch_versions = ch_versions.mix(CREATE_UMI_CONSENSUS.out.versions)
        } else {
            ch_reads_fastp = ch_input_fastq
        }

        // Trimming and/or splitting
        if (params.trim_fastq || params.split_fastq > 0) {

            save_trimmed_fail = false
            save_merged = false
            FASTP(ch_reads_fastp, save_trimmed_fail, save_merged)

            ch_reports = ch_reports.mix(
                                    FASTP.out.json.collect{meta, json -> json},
                                    FASTP.out.html.collect{meta, html -> html}
                                    )

            if(params.split_fastq){
                ch_reads_to_map = FASTP.out.reads.map{ key, reads ->

                        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                        [[
                            data_type:key.data_type,
                            id:key.id,
                            numLanes:key.numLanes,
                            patient: key.patient,
                            read_group:key.read_group,
                            sample:key.sample,
                            sex:key.sex,
                            size:read_files.size(),
                            status:key.status,
                        ],
                        read_files]
                    }.transpose()
            }else{
                ch_reads_to_map = FASTP.out.reads
            }

            ch_versions = ch_versions.mix(FASTP.out.versions)
        } else {
            ch_reads_to_map = ch_reads_fastp
        }

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        // reads will be sorted
        ch_reads_to_map = ch_reads_to_map.map{ meta, reads ->
            // update ID when no multiple lanes or splitted fastqs
            new_id = meta.size * meta.numLanes == 1 ? meta.sample : meta.id

            [[
                data_type:  meta.data_type,
                id:         new_id,
                numLanes:   meta.numLanes,
                patient:    meta.patient,
                read_group: meta.read_group,
                sample:     meta.sample,
                sex:        meta.sex,
                size:       meta.size,
                status:     meta.status,
                ],
            reads]
        }

        //
        // Logic to separate DNA from RNA samples, DNA samples will be aligned with bwa, and RNA samples with star
        //
        ch_reads_to_map.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
            rna:  it[0].status == 2
        }.set{ch_reads_to_map_status}





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

        // If no sex specified, sex is not considered
        // sex is only mandatory for somatic CNV
        // TODO remove as this is only specific to sarek, no CN is done in here
        if (row.sex) meta.sex = row.sex.toString()
        else meta.sex = 'NA'

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
            def tools_tumor = ['sage', 'controlfreec', 'mutect2', 'strelka2']  // This will be applied to tumour DNA and tumour RNA
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

            meta.numLanes   = numLanes.toInteger()
            meta.read_group = read_group.toString()
            meta.data_type  = 'bam'

            meta.size       = 1 // default number of splitted fastq

            if (params.step != 'annotate') return [meta, bam, bai]
            else {
                log.error "Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
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
                log.error "Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
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
                log.error "Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
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
                log.error "Samplesheet contains cram files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
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
                log.error "Samplesheet contains bam files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
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
                log.error "Samplesheet contains vcf files but step is `$params.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations"
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
