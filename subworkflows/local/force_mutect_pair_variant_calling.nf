//
// PAIRED VARIANT CALLING
//
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk4/tumor_normal_somatic_variant_calling/main'
include { RUN_MPILEUP as RUN_MPILEUP_NORMAL         } from '../nf-core/variantcalling/mpileup/main'
include { RUN_MPILEUP as RUN_MPILEUP_TUMOR          } from '../nf-core/variantcalling/mpileup/main'

workflow PAIR_VARIANT_CALLING_MUTECT2 {
    take:
        cram_pair                     // channel: [mandatory] cram
        dict                          // channel: [mandatory] dict
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        germline_resource             // channel: [optional]  germline_resource
        germline_resource_tbi         // channel: [optional]  germline_resource_tbi
        intervals                     // channel: [mandatory] intervals/target regions
        intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
        intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
        panel_of_normals              // channel: [optional]  panel_of_normals
        panel_of_normals_tbi          // channel: [optional]  panel_of_normals_tbi


    main:

    ch_versions          = Channel.empty()
    mutect2_vcf          = Channel.empty()

    // Remap channel with intervals
    cram_pair_intervals = cram_pair.combine(intervals)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, num_intervals ->
            // If no interval file provided (0) then add empty list
            intervals_new = num_intervals == 0 ? [] : intervals

            [[
                id:             meta.tumor_id + "_vs_" + meta.normal_id,
                normal_id:      meta.normal_id,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sex:            meta.sex,
                status:         meta.status,
                tumor_id:       meta.tumor_id,
                alleles:        meta.alleles
            ],
            normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_new]
        }


    cram_pair_mutect2 = cram_pair_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                            [meta, [normal_cram[0], tumor_cram[0]], [normal_crai, tumor_crai], intervals]
                        } // TODO: why are crams in a list?

    GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(
        cram_pair_mutect2,
        fasta,
        fasta_fai,
        dict,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi
    )

    mutect2_vcf = GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.filtered_vcf
    ch_versions = ch_versions.mix(GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.versions)



    emit:
    mutect2_vcf

    versions    = ch_versions
}
