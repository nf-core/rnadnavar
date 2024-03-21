/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { fromSamplesheet             } from 'plugin/nf-validation'
include { paramsSummaryMap            } from 'plugin/nf-validation'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_rnadnavar_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// Build the genome index and other reference files
include { SAMPLESHEET_TO_CHANNEL                                  } from '../subworkflows/local/samplesheet_to_channel/main'
include { PREPARE_REFERENCE_AND_INTERVALS                         } from '../subworkflows/local/prepare_reference_and_intervals/main'
include { PREPARE_INTERVALS as PREPARE_INTERVALS_FOR_REALIGNMENT  } from '../subworkflows/local/prepare_intervals/main'
// Download annotation cache if needed
include { ENSEMBLVEP_DOWNLOAD               } from '../modules/nf-core/ensemblvep/download/main'

// Alignment
include { BAM_ALIGN                         } from '../subworkflows/local/bam_align/main'

// Core subworkflows of the pipeline
include { BAM_VARIANT_CALLING_PRE_POST_PROCESSING as BAM_PROCESSING } from '../subworkflows/local/bam_variant_calling_pre_post_processing/main'

// Second run
include { BAM_EXTRACT_READS_HISAT2_ALIGN as PREPARE_REALIGNMENT  } from '../subworkflows/local/prepare_second_run/main'
include { BAM_VARIANT_CALLING_PRE_POST_PROCESSING as REALIGNMENT } from '../subworkflows/local/bam_variant_calling_pre_post_processing/main'

// Filter RNA
include { MAF_FILTERING_RNA } from '../subworkflows/local/maf_rna_filtering/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNADNAVAR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    // To gather all QC reports for MultiQC
    multiqc_files = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

    // Set input, can either be from --input or from automatic retrieval in utils_nfcore_rnadnavar_pipeline/main
    if (params.input) {
        ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet("input")
    } else {
        ch_from_samplesheet = params.build_only_index ? Channel.empty() : Channel.fromSamplesheet("input_restart")
    }
    // Parse samplesheet
    SAMPLESHEET_TO_CHANNEL(ch_from_samplesheet)

    input_sample = SAMPLESHEET_TO_CHANNEL.out.input_sample


    // Initialise MULTIQC
    multiqc_config                        = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    multiqc_custom_config                 = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    multiqc_logo                          = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

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

    multiqc_files = multiqc_files.mix(BAM_ALIGN.out.reports)
    versions = versions.mix(BAM_ALIGN.out.versions)

    // 5 MAIN STEPS: GATK PREPROCESING - VARIANT CALLING - NORMALIZATION - CONSENSUS - ANNOTATION
    BAM_PROCESSING(
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
        null,  // to repeat rescue consensus
        null,  // to repeat rescue consensus
        false,  // is second run
        params.no_intervals
    )
    filtered_maf = BAM_PROCESSING.out.maf
    multiqc_files      = multiqc_files.mix(BAM_PROCESSING.out.reports)
    versions     = versions.mix(BAM_PROCESSING.out.versions)
    if (params.tools && params.tools.split(',').contains('second_run')) {
        // fastq will not be split when realignment
        params.split_fastq = 0
        // reset intervals to none (realignment files are small)
        PREPARE_INTERVALS_FOR_REALIGNMENT(fasta_fai, null, true)

        PREPARE_REALIGNMENT(
                            input_sample,           // input from CSV if applicable
                            filtered_maf,
                            BAM_PROCESSING.out.cram_variant_calling,  // input from mapping
                            fasta,
                            fasta_fai,
                            dict,
                            PREPARE_REFERENCE_AND_INTERVALS.out.hisat2_index,
                            PREPARE_REFERENCE_AND_INTERVALS.out.splicesites,
                            BAM_PROCESSING.out.dna_consensus_maf,
                            BAM_PROCESSING.out.dna_varcall_mafs
                            ) // do mapping with hisat2

        versions = versions.mix(PREPARE_REALIGNMENT.out.versions)

        REALIGNMENT(
            Channel.empty(),                   // input from CSV if applicable: already processed in previous subworkflow
            PREPARE_REALIGNMENT.out.bam_mapped, // input from mapping
            Channel.empty(),                  // no cram from hisat2 for now
            fasta,
            fasta_fai,
            dict,
            dbsnp,
            dbsnp_tbi,
            pon,
            pon_tbi,
            known_sites_indels,
            known_sites_indels_tbi,
            germline_resource,
            germline_resource_tbi,
            PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed, // [[],0]
            PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_for_preprocessing, // [[id],0]
            PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_gz_tbi,  // [[[],[]],0]
            PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_combined,  // []
            PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_and_num_intervals,  // [[], 0]
            PREPARE_INTERVALS_FOR_REALIGNMENT.out.intervals_bed_gz_tbi_combined,  //[[],[]]
            PREPARE_REALIGNMENT.out.dna_consensus_maf,  // to repeat rescue consensus
            PREPARE_REALIGNMENT.out.dna_varcall_mafs,   // to repeat rescue consensus
            true,  // is realignment
            true // no_intervals
            )

        multiqc_files                = multiqc_files.mix(REALIGNMENT.out.reports)
        versions               = versions.mix(REALIGNMENT.out.versions)
        realigned_filtered_maf = REALIGNMENT.out.maf
    } else{
        realigned_filtered_maf = Channel.empty()
    }
    filtered_maf.dump(tag:"filtered_maf")
    realigned_filtered_maf.dump(tag:"realigned_filtered_maf")
    MAF_FILTERING_RNA(
                    filtered_maf.branch{  dna:  it[0].status < 2
                                            rna:  it[0].status >= 2
                                        }.rna,
                    realigned_filtered_maf.branch{  dna:  it[0].status < 2
                                            rna:  it[0].status >= 2
                                        }.rna,
                    fasta,
                    fasta_fai,
                    input_sample
                    )
    versions = versions.mix(MAF_FILTERING_RNA.out.versions)
//
// REPORTING
//
    version_yaml = Channel.empty()
    if (!(params.skip_tools && params.skip_tools.split(',').contains('versions'))) {
        CUSTOM_DUMPSOFTWAREVERSIONS(versions.unique().collectFile(name: 'collated_versions.yml'))
        version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }

    if (!(params.skip_tools && params.skip_tools.split(',').contains('multiqc'))) {
        workflow_summary = paramsSummaryMultiqc(paramsSummaryMap(workflow))
        workflow_summary = Channel.value(workflow_summary)

        methods_description = methodsDescriptionText(ch_multiqc_custom_methods_description)
        methods_description = Channel.value(methods_description)

        multiqc_files = Channel.empty()
        multiqc_files = multiqc_files.mix(version_yaml)
        multiqc_files = multiqc_files.mix(workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        multiqc_files = multiqc_files.mix(methods_description.collectFile(name: 'methods_description_mqc.yaml'))
        multiqc_files = multiqc_files.mix(multiqc_files.collect().ifEmpty([]))

    MULTIQC (
        multiqc_files.collect(),
        multiqc_config.toList(),
        multiqc_custom_config.toList(),
        multiqc_logo.toList()
    )
    }

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = versions                 // channel: [ path(versions.yml) ]
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
