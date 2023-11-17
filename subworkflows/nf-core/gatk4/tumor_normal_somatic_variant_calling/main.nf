//
// Run GATK mutect2 in tumor normal mode, getepileupsummaries, calculatecontamination, learnreadorientationmodel and filtermutectcalls
//

include { GATK4_MERGEVCFS                 as MERGE_MUTECT2               } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { GATK4_CALCULATECONTAMINATION    as CALCULATECONTAMINATION      } from '../../../../modules/nf-core/modules/gatk4/calculatecontamination/main'
include { GATK4_FILTERMUTECTCALLS         as FILTERMUTECTCALLS           } from '../../../../modules/nf-core/modules/gatk4/filtermutectcalls/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_NORMAL} from '../../../../modules/nf-core/modules/gatk4/gatherpileupsummaries/main'
include { GATK4_GATHERPILEUPSUMMARIES     as GATHERPILEUPSUMMARIES_TUMOR } from '../../../../modules/nf-core/modules/gatk4/gatherpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_NORMAL   } from '../../../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_GETPILEUPSUMMARIES        as GETPILEUPSUMMARIES_TUMOR    } from '../../../../modules/nf-core/modules/gatk4/getpileupsummaries/main'
include { GATK4_LEARNREADORIENTATIONMODEL as LEARNREADORIENTATIONMODEL   } from '../../../../modules/nf-core/modules/gatk4/learnreadorientationmodel/main'
include { GATK4_MERGEMUTECTSTATS          as MERGEMUTECTSTATS            } from '../../../../modules/nf-core/modules/gatk4/mergemutectstats/main'
include { GATK4_MUTECT2                   as MUTECT2                     } from '../../../../modules/nf-core/modules/gatk4/mutect2/main'

workflow GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING {
    take:
        input                     // channel: [ val(meta), [ input ], [ input_index ], [which_norm] ]
        fasta                     // channel: /path/to/reference/fasta
        fai                       // channel: /path/to/reference/fasta/index
        dict                      // channel: /path/to/reference/fasta/dictionary
        germline_resource         // channel: /path/to/germline/resource
        germline_resource_tbi     // channel: /path/to/germline/index
        panel_of_normals          // channel: /path/to/panel/of/normals
        panel_of_normals_tbi      // channel: /path/to/panel/of/normals/index
        skip_tools
        prior_contamination
        prior_segmentation
        prior_orientation

    main:
        ch_versions = Channel.empty()

        ch_input = input.map{ meta, input, input_index, intervals ->
                            [[
                            id: meta.tumor_id + "_vs_" + meta.normal_id,
                                            normal_id:      meta.normal_id,
                                            num_intervals:  meta.num_intervals,
                                            patient:        meta.patient,
                                            status:         meta.status,
                                            tumor_id:       meta.tumor_id,
                                            alleles:        meta.alleles
                            ], input, input_index, intervals]}
        MUTECT2(ch_input,
                fasta,
                fai,
                dict,
                germline_resource,
                germline_resource_tbi,
                panel_of_normals,
                panel_of_normals_tbi)
        // Figure out if using intervals or no_intervals
        MUTECT2.out.vcf.branch{
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }.set{ mutect2_vcf_branch }

        MUTECT2.out.tbi.branch{
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }.set{ mutect2_tbi_branch }

        MUTECT2.out.stats.branch{
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }.set{ mutect2_stats_branch }

        MUTECT2.out.f1r2.branch{
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }.set{ mutect2_f1r2_branch }
        mutect2_vcf_branch.intervals.dump(tag:'[STEP3] GATK4 - mutect2_vcf_branch.intervals')
        mutect2_to_merge = mutect2_vcf_branch.intervals
                                            .unique()  // Potential bug is nextflow were some files are duplicated
                                            .groupTuple()
        mutect2_to_merge.dump(tag:'[STEP3] GATK4 - mutect2_to_merge')
        //Only when using intervals
        MERGE_MUTECT2(
            mutect2_to_merge,
            dict
        )
        mutect2_vcf = Channel.empty().mix(
            MERGE_MUTECT2.out.vcf,
            mutect2_vcf_branch.no_intervals)
        mutect2_vcf.dump(tag:'[STEP3] GATK4 - mutect2')
        mutect2_tbi = Channel.empty().mix(
            MERGE_MUTECT2.out.tbi,
            mutect2_tbi_branch.no_intervals)
        mutect2_tbi.dump(tag:'[STEP3] GATK4 - mutect2 tbi')
        //Merge Mutect2 Stats
        MERGEMUTECTSTATS(
            mutect2_stats_branch.intervals
            .groupTuple())

        mutect2_stats = MERGEMUTECTSTATS.out.stats.mix(
            mutect2_stats_branch.no_intervals)
        mutect2_stats.dump(tag:'[STEP3] GATK4 - MUTECTSTATS')

        if (prior_orientation) {
            artifactprior = mutect2_vcf.map{meta, vcf -> [meta, file("${meta.id}_PRIOR_NO_ARTPRIOR")]}
        }
        else {

            if (!(skip_tools && skip_tools.split(',').contains('learnreadorientation'))){
            //Generate artifactpriors using learnreadorientationmodel on the f1r2 output of mutect2.
            LEARNREADORIENTATIONMODEL(mutect2_f1r2_branch.intervals
                                                        .groupTuple()
                                                        .mix(mutect2_f1r2_branch.no_intervals)
                                                        )
            artifactprior = LEARNREADORIENTATIONMODEL.out.artifactprior
            ch_versions   = ch_versions.mix(LEARNREADORIENTATIONMODEL.out.versions)
            }
            else{
                artifactprior = mutect2_vcf.map{meta, vcf -> [meta, file("${meta.id}_NO_ARTPRIOR")]}
            }
        }

        if (prior_contamination) {
            // prior contamination and segmentation tables calculated from a previous run
            segmentation = mutect2_vcf.map{meta, vcf -> [meta, file("${meta.id}_PRIOR_NO_SEG")]}
            contamination =  mutect2_vcf.map{meta, vcf -> [meta, file("${meta.id}_PRIOR_NO_TABLE")]}
            gather_table_tumor  = Channel.empty()
            gather_table_normal = Channel.empty()
        } else {
            if (!(skip_tools && skip_tools.split(',').contains('contamination'))){
                    //Generate pileup summary tables using getepileupsummaries. tumor sample should always be passed in as the first input and input list entries of ch_mutect2_in,
                    //to ensure correct file order for calculatecontamination.
                    pileup = input.multiMap{  meta, input_list, input_index_list, intervals ->
                        tumor: [ meta, input_list[1], input_index_list[1], intervals ]
                        normal: [ meta, input_list[0], input_index_list[0], intervals ]
                    }
                    germline_resource_pileup = germline_resource_tbi ? germline_resource : Channel.empty()
                    germline_resource_pileup_tbi = germline_resource_tbi ?: Channel.empty()
                    GETPILEUPSUMMARIES_TUMOR ( pileup.tumor.map{
                                                    meta, cram, crai, intervals ->
                                                    [[
                                                        id:             meta.tumor_id,
                                                        normal_id:      meta.normal_id,
                                                        num_intervals:  meta.num_intervals,
                                                        patient:        meta.patient,
                                                        status:         meta.status,
                                                        tumor_id:       meta.tumor_id,
                                                        alleles:        meta.alleles
                                                    ],
                                                        cram, crai, intervals]
                                                },
                                                fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi )
                    GETPILEUPSUMMARIES_NORMAL ( pileup.normal.map{
                                                    meta, cram, crai, intervals ->
                                                    [[
                                                        id:             meta.normal_id,
                                                        normal_id:      meta.normal_id,
                                                        num_intervals:  meta.num_intervals,
                                                        patient:        meta.patient,
                                                        status:         meta.status,
                                                        tumor_id:       meta.tumor_id,
                                                        alleles:        meta.alleles
                                                    ],
                                                        cram, crai, intervals]},
                                                fasta, fai, dict, germline_resource_pileup, germline_resource_pileup_tbi )
                    GETPILEUPSUMMARIES_NORMAL.out.table.branch{
                            intervals:    it[0].num_intervals > 1
                            no_intervals: it[0].num_intervals <= 1
                        }set{ pileup_table_normal }

                    GETPILEUPSUMMARIES_TUMOR.out.table.branch{
                            intervals:    it[0].num_intervals > 1
                            no_intervals: it[0].num_intervals <= 1
                        }set{ pileup_table_tumor }
                    //Merge Pileup Summaries
                    GATHERPILEUPSUMMARIES_NORMAL(
                        GETPILEUPSUMMARIES_NORMAL.out.table.groupTuple(),
                        dict)
                    gather_table_normal = Channel.empty().mix(
                        GATHERPILEUPSUMMARIES_NORMAL.out.table,
                        pileup_table_normal.no_intervals).map{ meta, table ->
                            new_meta = [
                                            id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                            normal_id:      meta.normal_id,
                                            num_intervals:  meta.num_intervals,
                                            patient:        meta.patient,
                                            status:         meta.status,
                                            tumor_id:       meta.tumor_id,
                                            alleles:        meta.alleles
                                        ]
                            [new_meta, table]
                        }
                    GATHERPILEUPSUMMARIES_TUMOR(
                        GETPILEUPSUMMARIES_TUMOR.out.table.groupTuple(),
                        dict
                        )
                    gather_table_tumor = Channel.empty().mix(
                        GATHERPILEUPSUMMARIES_TUMOR.out.table,
                        pileup_table_tumor.no_intervals).map{ meta, table ->
                            new_meta = [
                                        id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                        normal_id:      meta.normal_id,
                                        num_intervals:  meta.num_intervals,
                                        patient:        meta.patient,
                                        status:         meta.status,
                                        tumor_id:       meta.tumor_id,
                                        alleles:        meta.alleles
                                    ]
                            [new_meta, table]
                        }
                    //Contamination and segmentation tables created using calculatecontamination on the pileup summary table.
                    CALCULATECONTAMINATION( gather_table_tumor.join(gather_table_normal) )
                    segmentation = CALCULATECONTAMINATION.out.segmentation
                    contamination = CALCULATECONTAMINATION.out.contamination

                    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_NORMAL.out.versions)
                    ch_versions = ch_versions.mix(GETPILEUPSUMMARIES_TUMOR.out.versions)
                    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_NORMAL.out.versions)
                    ch_versions = ch_versions.mix(GATHERPILEUPSUMMARIES_TUMOR.out.versions)
                    ch_versions = ch_versions.mix(CALCULATECONTAMINATION.out.versions)
            }
            else {
                segmentation = mutect2_vcf.map{meta, vcf -> [meta, file("${meta.id}_NO_SEG")]}
                contamination =  mutect2_vcf.map{meta, vcf -> [meta, file("${meta.id}_NO_TABLE")]}
                gather_table_tumor = Channel.empty()
                gather_table_normal = Channel.empty()
            }
        }
        mutect2_vcf.dump(tag:"mutect2_vcf FOR FMC")
        contamination.dump(tag:"CONTAMINATION FOR FMC")
        artifactprior.dump(tag:"artifactprior FOR FMC")
        segmentation.dump(tag:"segmentation FOR FMC")
        ch_filtermutect    = mutect2_vcf.join(mutect2_tbi)
                                        .join(mutect2_stats)
                                        .join(artifactprior)
                                        .join(segmentation)
                                        .join(contamination)

        //Mutect2 calls filtered by filtermutectcalls using the artifactpriors, contamination and segmentation tables.
        ch_filtermutect_in = ch_filtermutect.map{ meta, vcf, tbi, stats, orientation, seg, cont ->
                                                    [meta, vcf, tbi, stats, orientation, seg, cont, []] }
        ch_filtermutect.dump(tag:"[STEP3] GATK4 - input FILTERMUTECTCALLS")
        FILTERMUTECTCALLS ( ch_filtermutect_in, fasta, fai, dict )
        filtered_vcf = FILTERMUTECTCALLS.out.vcf.map{ meta, vcf ->
                                                [[patient:meta.patient,
                                                normal_id:meta.normal_id,
                                                tumor_id:meta.tumor_id,
                                                status:meta.status,
                                                id:meta.tumor_id + "_vs_" + meta.normal_id,
                                                num_intervals:meta.num_intervals,
                                                alleles:meta.alleles,
                                                variantcaller:"mutect2"
                                                ], vcf]}
        ch_versions = ch_versions.mix(MERGE_MUTECT2.out.versions)
        ch_versions = ch_versions.mix(FILTERMUTECTCALLS.out.versions)


        ch_versions = ch_versions.mix(MERGEMUTECTSTATS.out.versions)
        ch_versions = ch_versions.mix(MUTECT2.out.versions)

    emit:
        mutect2_vcf            = mutect2_vcf                                    // channel: [ val(meta), [ vcf ] ]
        mutect2_stats          = mutect2_stats                                  // channel: [ val(meta), [ stats ] ]
        artifact_priors        = artifactprior   // channel: [ val(meta), [ artifactprior ] ]
        pileup_table_tumor     = gather_table_tumor                             // channel: [ val(meta), [ table_tumor ] ]
        pileup_table_normal    = gather_table_normal                            // channel: [ val(meta), [ table_normal ] ]
        contamination_table    = contamination       // channel: [ val(meta), [ contamination ] ]
        segmentation_table     = segmentation        // channel: [ val(meta), [ segmentation ] ]
        filtered_vcf           = filtered_vcf                                   // channel: [ val(meta), [ vcf ] ]
        filtered_tbi           = FILTERMUTECTCALLS.out.tbi                      // channel: [ val(meta), [ tbi ] ]
        filtered_stats         = FILTERMUTECTCALLS.out.stats                    // channel: [ val(meta), [ stats ] ]

        versions               = ch_versions                                    // channel: [ versions.yml ]
}
