//
// STEP7 : FILTERING VARIANTS
//

include { BASIC_FILTERING                               } from '../../../modules/local/filter_variants'
include { BASIC_FILTERING as BASIC_FILTERING_RNA        } from '../../../modules/local/filter_variants'
include { SECOND_PASS as SECOND_PASS_RNA                } from './second_pass'
include { VCF2MAF                                       } from '../../../modules/local/vcf2maf/vcf2maf/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_SECOND_PASS  } from '../../../modules/nf-core/modules/samtools/merge/main'
include { RNA_FILTERING                                 } from '../../../modules/local/rna_filtering'



workflow FILTERING {

    take:
        ch_input_sample
        vcf_to_maf
        vcf_consensus_dna // to repeat rescue consensus
        vcfs_status_dna // to repeat rescue consensus
        star_bams
        bwa_bams
        fasta
        fasta_fai
        dict
        hisat2_index
        splicesites
        dbsnp
        dbsnp_tbi
        pon
        pon_tbi
        germline_resource
        germline_resource_tbi
        intervals
        intervals_for_preprocessing
        ch_interval_list_split
        intervals_bed_gz_tbi
        intervals_bed_combined


    main:
        ch_reports   = Channel.empty()
        ch_versions  = Channel.empty()

        if (params.step == 'filtering')  vcf_to_maf = ch_input_sample

        // Initiate resources variables
        whitelist  = params.whitelist  ? Channel.fromPath(params.whitelist).collect() : Channel.value([])
        darned     = params.darned     ? Channel.fromPath(params.darned).collect()    : Channel.value([])
        radar      = params.radar      ? Channel.fromPath(params.radar).collect()     : Channel.value([])
        nat        = params.nat        ? Channel.fromPath(params.nat).collect()       : Channel.value([])
        redi       = params.redi       ? Channel.fromPath(params.redi).collect()      : Channel.value([])


        // First we transform the vcf to MAF
        VCF2MAF(vcf_to_maf,
                fasta)
        maf_to_filter = VCF2MAF.out.maf
        maf_to_filter.dump(tag:"[STEP7: FILTERING] maf input")
        whitelist.dump(tag:"[STEP7: FILTERING] whitelist")
        ch_versions = ch_versions.mix(VCF2MAF.out.versions)

        // BASIC FILTERING
        BASIC_FILTERING(maf_to_filter, whitelist, fasta)
        BASIC_FILTERING.out.maf.dump(tag:"[STEP7: FILTERING] maf filtered")
        // Once this is done DNA is ready, RNA still has a 2nd PASS TODO: optional?
        BASIC_FILTERING.out.maf.branch{
                dna: it[0].status < 2
                rna: it[0].status == 2
            }.set{maf_to_filter_status}
        ch_versions = ch_versions.mix(BASIC_FILTERING.out.versions)

        // STEP 8: RNA FILTERING - TODO: make it optional and run just rna filtering
        // RNA specific filtering (2nd PASS) - this fast BUT it increases the length of the pipeline considerably

        // 1 We take the previous aligned reads with star for tumor RNA and DNA normal
        bwa_bams.branch{
                        normal: it[0].status == 0
                        tumor: it[0].status == 1
                    }.set{previous_dna_alignment}
        // we only need normals - dna tumour will NOT be realigned
        previous_normal_alignment = previous_dna_alignment.normal.groupTuple()
        // 2. Group them and merge if applicable
        previous_alignment = star_bams
                                    .mix(previous_normal_alignment)
        SAMTOOLS_MERGE_SECOND_PASS(previous_alignment, fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SECOND_PASS.out.versions)

        previous_alignment_merged = SAMTOOLS_MERGE_SECOND_PASS.out.bam
                                        .map{meta, bam -> [
                                                            [
                                                                id:meta.sample,
                                                                data_type:"bam",
                                                                patient:meta.patient,
                                                                sample:meta.sample,
                                                                read_group:meta.read_group,
                                                                status:meta.status
                                                            ],
                                                            bam
                                                          ]
                                             }
        previous_alignment_merged.dump(tag:"[STEP7: FILTERING] bams for realignment")
        SECOND_PASS_RNA(
                ch_input_sample,
                maf_to_filter_status.rna,
                previous_alignment_merged,
                vcf_consensus_dna,
                vcfs_status_dna,
                fasta,
                fasta_fai,
                dict,
                hisat2_index,
                splicesites,
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
                intervals_bed_combined
                )
        ch_versions = ch_versions.mix(SECOND_PASS_RNA.out.versions)
        ch_reports = ch_versions.mix(SECOND_PASS_RNA.out.reports)

        maf_to_filter = SECOND_PASS_RNA.out.maf // TODO: optional?
        maf_to_filter.dump(tag:"[STEP7: FILTERING] maf_to_filter_rna_2pass")
        BASIC_FILTERING_RNA(maf_to_filter, whitelist, fasta)

        maf_to_filter_status.rna.dump(tag:"[STEP7: FILTERING] maf_to_filter_status.rna")
        BASIC_FILTERING_RNA.out.maf.dump(tag:"[STEP7: FILTERING] BASIC_FILTERING_RNA.out.maf")

        maf_to_cross_first_pass = maf_to_filter_status.rna
                                    .map{meta, maf -> [meta.patient, meta, maf]}
        maf_to_cross_second_pass = BASIC_FILTERING_RNA.out.maf
                                    .map{meta, maf -> [meta.patient, meta, maf]}

        maf_to_cross_first_pass.dump(tag:"[STEP7: FILTERING] maf_to_cross_first_pass")
        maf_to_cross_second_pass.dump(tag:"[STEP7: FILTERING] maf_to_cross_second_pass")
        maf_to_cross_first_pass
                        .cross(maf_to_cross_second_pass).dump(tag:"[STEP7: FILTERING] maf_to_cross_crossed_pass")

        maf_crossed = maf_to_cross_first_pass
                        .cross(maf_to_cross_second_pass)
                        .map{first, second ->
                                def meta = [:]
                                meta.patient    = first[0]
                                meta.first_id   = first[1].tumor_id
                                meta.second_id  = second[1].tumor_id
                                meta.status     = first[1].status
                                meta.tumor_id   = first[1].tumor_id
                                meta.id   = first[1].tumor_id
                                meta.normal_id  = first[1].normal_id
                                [meta, first[2], second[2]]
                                }
        maf_crossed.dump(tag:"[STEP7: FILTERING] maf_crossed")
//        maf_to_filter_status_dna = maf_to_filter_status.dna.map{meta, maf -> [meta, maf, file("NO_FILE.maf")]}
//        maf_to_filter_status_dna.dump(tag:"[STEP7: FILTERING] maf_to_filter_status_dna")
//        maf_crossed = maf_crossed.mix(maf_to_filter_status.dna)
        RNA_FILTERING(maf_crossed,
                      fasta)
        ch_versions = ch_versions.mix(RNA_FILTERING.out.versions)

        // TODO RNA PON
        // TODO produce some stats DNA vs RNA, oncoprint, etc
//        FINAL_REPORT(RNA_FILTERING.out.maf,
//                     BASIC_FILTERING.out,
//                     CONSENSUS.out.txt)





    emit:
        versions                         = ch_versions                                                         // channel: [ versions.yml ]
        reports                         = ch_reports                                                         // channel: [ versions.yml ]
}