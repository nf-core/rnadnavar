//
// PAIRED VARIANT CALLING
//
include { GATK_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING } from '../../subworkflows/nf-core/gatk4/tumor_normal_somatic_variant_calling/main'
include { RUN_MANTA_SOMATIC                         } from '../nf-core/variantcalling/manta/somatic/main.nf'
include { RUN_FREEBAYES as RUN_FREEBAYES_SOMATIC    } from '../nf-core/variantcalling/freebayes/main.nf'
include { RUN_STRELKA_SOMATIC                       } from '../nf-core/variantcalling/strelka/somatic/main.nf'
include { RUN_MPILEUP as RUN_MPILEUP_NORMAL         } from '../nf-core/variantcalling/mpileup/main'
include { RUN_MPILEUP as RUN_MPILEUP_TUMOR          } from '../nf-core/variantcalling/mpileup/main'

workflow PAIR_VARIANT_CALLING {
    take:
        tools                         // Mandatory, list of tools to apply
        cram_pair                     // channel: [mandatory] cram
        dbsnp                         // channel: [mandatory] dbsnp
        dbsnp_tbi                     // channel: [mandatory] dbsnp_tbi
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

    //TODO: Temporary until the if's can be removed and printing to terminal is prevented with "when" in the modules.config
    manta_vcf            = Channel.empty()
    strelka_vcf          = Channel.empty()
    mutect2_vcf          = Channel.empty()
    freebayes_vcf          = Channel.empty()

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
                tumor_id:       meta.tumor_id
            ],
            normal_cram, normal_crai, tumor_cram, tumor_crai, intervals_new]
        }
    // Remap channel with gzipped intervals + indexes
    cram_pair_intervals_gz_tbi = cram_pair.combine(intervals_bed_gz_tbi)
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed_tbi, num_intervals ->

            //If no interval file provided (0) then add empty list
            bed_new = num_intervals == 0 ? [] : bed_tbi[0]
            tbi_new = num_intervals == 0 ? [] : bed_tbi[1]

            [[
                id:             meta.tumor_id + "_vs_" + meta.normal_id,
                normal_id:      meta.normal_id,
                num_intervals:  num_intervals,
                patient:        meta.patient,
                sex:            meta.sex,
                status:         meta.status,
                tumor_id:       meta.tumor_id
            ],
            normal_cram, normal_crai, tumor_cram, tumor_crai, bed_new, tbi_new]

        }

    if (tools.split(',').contains('manta')) {
        RUN_MANTA_SOMATIC(
            cram_pair_intervals_gz_tbi,
            dict,
            fasta,
            fasta_fai
        )
        manta_vcf                            = RUN_MANTA_SOMATIC.out.manta_vcf
        manta_candidate_small_indels_vcf     = RUN_MANTA_SOMATIC.out.manta_candidate_small_indels_vcf
        manta_candidate_small_indels_vcf_tbi = RUN_MANTA_SOMATIC.out.manta_candidate_small_indels_vcf_tbi
        ch_versions                          = ch_versions.mix(RUN_MANTA_SOMATIC.out.versions)
    }

    if (tools.split(',').contains('strelka')) {

        if (tools.split(',').contains('manta')) {
            cram_pair_strelka = cram_pair.join(manta_candidate_small_indels_vcf)
                                         .join(manta_candidate_small_indels_vcf_tbi)
                                         .combine(intervals_bed_gz_tbi)
                                         .map{
                                            meta, normal_cram, normal_crai, tumor_cram, tumor_crai, vcf, vcf_tbi, bed_tbi, num_intervals ->
                                            //If no interval file provided (0) then add empty list
                                            bed_new = num_intervals <= 1 ? [] : bed_tbi[0]
                                            tbi_new = num_intervals <= 1 ? [] : bed_tbi[1]

                                            [
                                                [
                                                    id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                                    normal_id:      meta.normal_id,
                                                    num_intervals:  num_intervals,
                                                    patient:        meta.patient,
                                                    sex:            meta.sex,
                                                    status:         meta.status,
                                                    tumor_id:       meta.tumor_id
                                                ],
                                                normal_cram, normal_crai, tumor_cram, tumor_crai, vcf, vcf_tbi, bed_new, tbi_new
                                            ]
                                         }
        } else {
            cram_pair_strelka = cram_pair_intervals_gz_tbi.map{
                    meta, normal_cram, normal_crai, tumor_cram, tumor_crai, bed, tbi ->
                    [meta, normal_cram, normal_crai, tumor_cram, tumor_crai, [], [], bed, tbi]
            }
        }

        RUN_STRELKA_SOMATIC(
            cram_pair_strelka,
            dict,
            fasta,
            fasta_fai
            )

        strelka_vcf = Channel.empty().mix(RUN_STRELKA_SOMATIC.out.strelka_vcf)
        ch_versions = ch_versions.mix(RUN_STRELKA_SOMATIC.out.versions)
    }

    if (tools.split(',').contains('freebayes')){

        RUN_FREEBAYES_SOMATIC(
            cram_pair_intervals,
            dict,
            fasta,
            fasta_fai
        )

        freebayes_vcf = RUN_FREEBAYES_SOMATIC.out.freebayes_vcf
        ch_versions   = ch_versions.mix(RUN_FREEBAYES_SOMATIC.out.versions)
    }

    if (tools.split(',').contains('mutect2')) {
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
    }


    emit:
    manta_vcf
    strelka_vcf
    freebayes_vcf
    mutect2_vcf

    versions    = ch_versions
}
