//
// Filtering steps specific for RNA
//
include { RNA_FILTERING                                 } from '../../modules/local/rna_filtering'

workflow FILTERING_RNA {
    take:
        tools
        maf_to_filter_rna_1   // maf from first pass
        maf_to_filter_rna_2   // maf from realignment (second pass) [OPT?]
        fasta

    main:
        ch_versions = Channel.empty()
        if (tools.split(',').contains('rna_filtering')) {
            // RNA filtering after realignment
            maf_to_cross_first_pass = maf_to_filter_rna_1
                                        .map{meta, maf -> [meta.patient, meta, maf]}


            maf_to_cross_first_pass.dump(tag:"[STEP9: FILTERING] maf_to_cross_first_pass")

            if (maf_to_filter_rna_2){
                maf_to_cross_second_pass = maf_to_filter_rna_2
                                            .map{meta, maf -> [meta.patient, meta, maf]}
                maf_to_cross_second_pass.dump(tag:"[STEP9: FILTERING] maf_to_cross_second_pass")
                maf_to_cross_first_pass
                                .cross(maf_to_cross_second_pass).dump(tag:"[STEP9: FILTERING] maf_to_cross_crossed_pass")
                maf_crossed = maf_to_cross_first_pass
                                .cross(maf_to_cross_second_pass)
                                .map{first, second ->
                                        def meta = [:]
                                        meta.patient    = first[0]
                                        meta.first_id   = first[1].id
                                        meta.second_id  = second[1].id
                                        meta.status     = first[1].status
                                        meta.tumor_id   = first[1].rna_id.split("_vs_")[0]
                                        meta.id         = first[1].id
                                        meta.normal_id  = first[1].rna_id.split("_vs_")[1]
                                        [meta, first[2], second[2]]
                                        }
            } else {
            maf_crossed = maf_to_cross_first_pass.map{meta, maf -> [meta, maf, []]}

            }
            maf_crossed.dump(tag:"[STEP9: FILTERING] maf_crossed")
    //        maf_to_filter_status_dna = maf_to_filter_status.dna.map{meta, maf -> [meta, maf, file("NO_FILE.maf")]}
    //        maf_to_filter_status_dna.dump(tag:"[STEP9: FILTERING] maf_to_filter_status_dna")
    //        maf_crossed = maf_crossed.mix(maf_to_filter_status.dna)
            RNA_FILTERING(maf_crossed,
                          fasta)
            ch_versions = ch_versions.mix(RNA_FILTERING.out.versions)
        }


    emit:
        versions            = ch_versions // channel: [ versions.yml ]



}
