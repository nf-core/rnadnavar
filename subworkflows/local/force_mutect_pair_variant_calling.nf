//
// PAIRED VARIANT CALLING
//
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk4/tumor_normal_somatic_variant_calling/main'
include { RUN_MPILEUP as RUN_MPILEUP_NORMAL         } from '../nf-core/variantcalling/mpileup/main'
include { RUN_MPILEUP as RUN_MPILEUP_TUMOR          } from '../nf-core/variantcalling/mpileup/main'
include { VT_DECOMPOSE                              } from '../../modules/local/vt/decompose/main'
include { VCF_QC                                    } from '../nf-core/vcf_qc'


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
        skip_tools
        contamination
        segmentation
        orientation


    main:

        ch_versions          = Channel.empty()
        ch_reports           = Channel.empty()
        mutect2_vcf          = Channel.empty()

        // Remap channel with intervals
        cram_pair.dump(tag:"FORCE_CALLS cram_pair")
        cram_pair_intervals = cram_pair.combine(intervals)
            .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals, num_intervals ->
                // If no interval file provided (0) then add empty list
                intervals_new = num_intervals == 0 ? [] : intervals
                [[
                    id:             meta.tumor_id + "_vs_" + meta.normal_id,
                    normal_id:      meta.normal_id,
                    num_intervals:  num_intervals,
                    patient:        meta.patient,
                    status:         meta.status,
                    tumor_id:       meta.tumor_id,
                    alleles:        meta.alleles
                ],
                normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_new]
            }
        if (contamination){
            contamination = contamination.combine(intervals)
            .map {meta, table, intervals, num_intervals ->
                  // If no interval file provided (0) then add empty list
                intervals_new = num_intervals == 0 ? [] : intervals
                [meta + [num_intervals:  num_intervals], table]
                  }
            segmentation = segmentation.combine(intervals)
            .map {meta, table, intervals, num_intervals ->
                  // If no interval file provided (0) then add empty list
                intervals_new = num_intervals == 0 ? [] : intervals
                [meta + [num_intervals:  num_intervals], table]
                  }
        }
        orientation = orientation.combine(intervals)
            .map {meta, table, intervals, num_intervals ->
                  // If no interval file provided (0) then add empty list
                intervals_new = num_intervals == 0 ? [] : intervals
                [meta + [num_intervals:  num_intervals], table]
                  }

        cram_pair_intervals.dump(tag:"FORCE_CALLS cram_pair_intervals")
        cram_pair_mutect2 = cram_pair_intervals.map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervals ->
                                    if (meta.num_intervals == 1){
                                        [meta, [normal_cram[0], tumor_cram[0]], [normal_crai, tumor_crai], intervals]
                                    } else{
                                        [meta, [normal_cram, tumor_cram], [normal_crai, tumor_crai], intervals]}
                                }

        cram_pair_mutect2.dump(tag:"FORCE_CALLS cram_pair_mutect2")
        cram_pair_mutect2_test = cram_pair_mutect2.join(contamination).join(segmentation).join(orientation)
                          .map{meta, pairs, pairs_idx, intervals, cont_table, seg_table, artprior ->
                                [meta+[cont:cont_table, seg:seg_table, orient:artprior], pairs, pairs_idx, intervals]}
        cram_pair_mutect2_test.dump(tag:"FORCE_CALLS cram_pair_mutect2_test")

        GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING(
            cram_pair_mutect2_test,
            fasta,
            fasta_fai,
            dict,
            germline_resource,
            germline_resource_tbi,
            panel_of_normals,
            panel_of_normals_tbi,
            skip_tools,
            true,
            true,
            true
        )
        forced_vcf = GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.filtered_vcf
        forced_vcf.dump(tag:"FORCE_CALLS forced_vcf")
        ch_versions = ch_versions.mix(GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING.out.versions)
        // Decompose to facilitate filtering
        VT_DECOMPOSE(forced_vcf)
        VT_DECOMPOSE.out.vcf.dump(tag:"FORCE_CALLS VT_DECOMPOSE.out.vcf")
        ch_versions = ch_versions.mix(VT_DECOMPOSE.out.versions)

        // Get QC
        VCF_QC(forced_vcf, intervals_bed_combined)
        ch_reports  = ch_reports.mix(VCF_QC.out.bcftools_stats.collect{meta, stats -> stats})
        ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_counts.collect{ meta, counts -> counts})
        ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_qual.collect{ meta, qual -> qual })
        ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_filter_summary.collect{meta, summary -> summary})

    emit:
        vcf         = VT_DECOMPOSE.out.vcf
        versions    = ch_versions
        reports     = ch_reports
}
