//
// PAIRED VARIANT CALLING
//
include { BAM_VARIANT_CALLING_SOMATIC_MANTA             } from '../bam_variant_calling_somatic_manta/main'
include { BAM_VARIANT_CALLING_SOMATIC_MUTECT2           } from '../bam_variant_calling_somatic_mutect2/main'
include { BAM_VARIANT_CALLING_SOMATIC_STRELKA           } from '../bam_variant_calling_somatic_strelka/main'
include { BAM_VARIANT_CALLING_SOMATIC_SAGE              } from '../bam_variant_calling_somatic_sage/main'

workflow BAM_VARIANT_CALLING_SOMATIC {
    take:
    tools                         // Mandatory, list of tools to apply
    cram                          // channel: [mandatory] cram
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    dict                          // channel: [mandatory] dict
    germline_resource             // channel: [optional]  germline_resource
    germline_resource_tbi         // channel: [optional]  germline_resource_tbi
    intervals                     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals
    intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
    intervals_bed_gz_tbi_combined // [] if no_intervals, else interval_bed_combined_gz, interval_bed_combined_gz_tbi
    panel_of_normals              // channel: [optional]  panel_of_normals
    panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi
    joint_mutect2                 // boolean: [mandatory] [default: false] run mutect2 in joint mode
	second_run

    main:
    versions          = Channel.empty()

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config    vcf_manta         = Channel.empty()
    vcf_manta         = Channel.empty()
    vcf_strelka       = Channel.empty()
    vcf_mutect2       = Channel.empty()
    vcf_sage          = Channel.empty()
    // SAGE
    if (tools.split(',').contains('sage') || second_run) {

        BAM_VARIANT_CALLING_SOMATIC_SAGE(
            cram,
            // Remap channel to match module/subworkflow
            dict,
            // Remap channel to match module/subworkflow
            fasta.map{ it -> [ [ id:'fasta' ], it ] },
            // Remap channel to match module/subworkflow
            fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] },
            intervals
        )

        vcf_sage   = BAM_VARIANT_CALLING_SOMATIC_SAGE.out.vcf
        versions   = versions.mix(BAM_VARIANT_CALLING_SOMATIC_SAGE.out.versions)
    }

    // MANTA
    if (tools.split(',').contains('manta')) {
        BAM_VARIANT_CALLING_SOMATIC_MANTA(
            cram,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi_combined
        )

        vcf_manta = BAM_VARIANT_CALLING_SOMATIC_MANTA.out.vcf
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.versions)
    }

    // STRELKA
    if (tools.split(',').contains('strelka') || second_run) {
        // Remap channel to match module/subworkflow
        cram_strelka = (tools.split(',').contains('manta')) ?
            cram.join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf, failOnDuplicate: true, failOnMismatch: true).join(BAM_VARIANT_CALLING_SOMATIC_MANTA.out.candidate_small_indels_vcf_tbi, failOnDuplicate: true, failOnMismatch: true) :
            cram.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [] ] }

        BAM_VARIANT_CALLING_SOMATIC_STRELKA(
            cram_strelka,
            // Remap channel to match module/subworkflow
            dict,
            fasta,
            fasta_fai,
            intervals_bed_gz_tbi
        )

        vcf_strelka = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.vcf)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_STRELKA.out.versions)
    }

    // MUTECT2
    if (tools.split(',').contains('mutect2') || second_run) {
        BAM_VARIANT_CALLING_SOMATIC_MUTECT2(
            // Remap channel to match module/subworkflow
            cram.map { meta, normal_cram, normal_crai, tumor_cram, tumor_crai -> [ meta, [ normal_cram, tumor_cram ], [ normal_crai, tumor_crai ] ] },
            // Remap channel to match module/subworkflow
            fasta.map{ it -> [ [ id:'fasta' ], it ] },
            // Remap channel to match module/subworkflow
            fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] },
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            intervals,
            joint_mutect2,
            second_run
        )

        vcf_mutect2                 = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.vcf_filtered
        contamination_table_mutect2 = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.contamination_table
        segmentation_table_mutect2  = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.segmentation_table
        artifact_priors_mutect2     = BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.artifact_priors
        versions                    = versions.mix(BAM_VARIANT_CALLING_SOMATIC_MUTECT2.out.versions)
    } else {

        contamination_table_mutect2 = Channel.empty()
		segmentation_table_mutect2  = Channel.empty()
		artifact_priors_mutect2     = Channel.empty()


    }

    vcf_all = Channel.empty().mix(
        vcf_mutect2,
        vcf_strelka,
        vcf_sage
    )

    emit:
    vcf_all
    vcf_manta
    vcf_mutect2
    vcf_strelka
    vcf_sage
    contamination_table_mutect2 = contamination_table_mutect2
	segmentation_table_mutect2  = segmentation_table_mutect2
	artifact_priors_mutect2     = artifact_priors_mutect2
    versions
}