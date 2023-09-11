//
// Core workflow of the RNA/DNA variant calling pipeline
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
include { BAM_GATK_PREPROCESSING                          } from '../bam_gatk_preprocessing/main'
// For now only matched supported
include { BAM_VARIANT_CALLING                             } from '../bam_variant_calling/main'
// Normalise VCFs
include { VCF_NORMALISE                                   } from '../vcf_normalise/main'
// Annotation
include { VCF_ANNOTATE                                    } from '../vcf_annotate/main'
// Consensus
 include { VCF_CONSENSUS                                  } from '../vcf_consensus/main'
// Filtering
 include { MAF_FILTERING                                  } from '../maf_filtering/main'


workflow BAM_VARIANT_CALLING_PRE_POST_PROCESSING {
    take:
    input_sample                    // input from CSV if applicable
    bam_mapped                      // channel: [mandatory] bam_mapped
    cram_mapped                     // channel: [mandatory] cram_mapped
    fasta                           // fasta reference file
    fasta_fai                       // fai for fasta file
    dict                            // dict for fasta file
    dbsnp                           // channel: [optional]  germline_resource
    dbsnp_tbi                       // channel: [optional]  germline_resource_tbi
    pon                             // channel: [optional]  pon for mutect2
    pon_tbi                         // channel: [optional]  pon_tbi for mutect2
    known_sites_indels              // channel: [optional]  known_sites
    known_sites_indels_tbi          // channel: [optional]  known_sites
    germline_resource               // channel: [optional]  germline_resource
    germline_resource_tbi           // channel: [optional]  germline_resource
    intervals                       // channel: [mandatory] intervals/target regions
    intervals_for_preprocessing     // channel: [mandatory] intervals/wes
    intervals_bed_gz_tbi            // channel: [mandatory] intervals/target regions index zipped and indexed
    intervals_bed_combined          // channel: [mandatory] intervals/target regions in one file unzipped
    intervals_and_num_intervals     // channel: [mandatory] [ intervals, num_intervals ] (or [ [], 0 ] if no intervals)
    intervals_bed_gz_tbi_combined   // channel: [mandatory] intervals/target regions in one file zipped
    dna_consensus_maf               // to repeat rescue consensus
    dna_varcall_mafs                // to repeat rescue consensus
    second_run

    main:
    reports   = Channel.empty()
    versions  = Channel.empty()

	// GATK PREPROCESSING - See: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
    BAM_GATK_PREPROCESSING(
        input_sample,
        bam_mapped,                               // channel: [mandatory] [meta, [bam]]
        cram_mapped,                              // channel: [mandatory] [meta, [cram]]
        fasta,                                    // channel: [mandatory] fasta
        fasta_fai ,                               // channel: [mandatory] fasta_fai
        dict.map{ it -> [ [ id:'dict' ], it[0] ] },  // channel: [mandatory] dict
	    known_sites_indels,                       // channel: [optional]  known_sites
	    known_sites_indels_tbi,                   // channel: [optional]  known_sites
        germline_resource,                        // channel: [optional]  germline_resource
        germline_resource_tbi,                    // channel: [optional]  germline_resource_tbi
        intervals,                                // channel: [mandatory] intervals/target regions
        intervals_for_preprocessing,              // channel: [mandatory] intervals_for_preprocessing/wes
        intervals_and_num_intervals,              // channel: [mandatory] intervals_for_preprocessing/wes
        second_run
    )

    cram_variant_calling = BAM_GATK_PREPROCESSING.out.cram_variant_calling
    versions             = versions.mix(BAM_GATK_PREPROCESSING.out.versions)
    reports              = reports.mix(BAM_GATK_PREPROCESSING.out.reports)

	// VARIANT CALLING
     BAM_VARIANT_CALLING(
         params.tools,
         cram_variant_calling,
         fasta,
         fasta_fai,
         dict,
         germline_resource,
         germline_resource_tbi,
         intervals,
         intervals_bed_gz_tbi,
         intervals_bed_combined,
         intervals_bed_gz_tbi_combined,
         pon,
         pon_tbi,
         input_sample,
         second_run
     )
     cram_variant_calling_pair     = BAM_VARIANT_CALLING.out.cram_variant_calling_pair  // use same crams for force calling later
     vcf_to_normalise              = BAM_VARIANT_CALLING.out.vcf_to_normalise
     contamination                 = BAM_VARIANT_CALLING.out.contamination_table
     segmentation                  = BAM_VARIANT_CALLING.out.segmentation_table
     orientation                   = BAM_VARIANT_CALLING.out.artifact_priors
     versions                      = versions.mix(BAM_VARIANT_CALLING.out.versions)
     reports                       = reports.mix(BAM_VARIANT_CALLING.out.reports)


    // NORMALISE
    VCF_NORMALISE (
                   vcf_to_normalise,
                   // Remap channel to match module/subworkflow
                   fasta.map{ it -> [ [ id:'fasta' ], it ] },
                   input_sample,
                   second_run
                   )
    versions                      = versions.mix(VCF_NORMALISE.out.versions)
    vcf_to_annotate               = VCF_NORMALISE.out.vcf

    // ANNOTATION

    VCF_ANNOTATE(
                vcf_to_annotate.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
                fasta,
                input_sample,
                second_run
                )

    vcf_to_consensus              = VCF_ANNOTATE.out.vcf_ann
    versions                      = versions.mix(VCF_ANNOTATE.out.versions)
    reports                       = reports.mix(VCF_ANNOTATE.out.reports)

	vcf_to_consensus.dump(tag:"vcf_to_consensus0")
	// STEP 6: CONSENSUS
	VCF_CONSENSUS (
			params.tools,
			vcf_to_consensus,
			fasta,
			dna_consensus_maf, // null when first pass
			dna_varcall_mafs,  // null when first pass
			input_sample,
            second_run
	         )

    dna_consensus_maf             = VCF_CONSENSUS.out.maf_consensus_dna
    dna_varcall_mafs              = VCF_CONSENSUS.out.mafs_dna
    maf_to_filter                 = VCF_CONSENSUS.out.maf
    versions                      = versions.mix(VCF_CONSENSUS.out.versions)

	maf_to_filter.dump(tag:"maf_to_filter0")
    // STEP 7: FILTERING
	MAF_FILTERING(maf_to_filter, fasta, input_sample, second_run)
	filtered_maf = MAF_FILTERING.out.maf
    versions     = versions.mix(MAF_FILTERING.out.versions)


     emit:
     dna_consensus_maf           = dna_consensus_maf
     dna_varcall_mafs            = dna_varcall_mafs
     maf                         = filtered_maf
     versions                    = versions  // channel: [ versions.yml ]
     reports                     = reports
}
