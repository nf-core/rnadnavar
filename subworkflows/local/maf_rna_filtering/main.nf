//
// Filtering steps specific for RNA
//
include { RNA_FILTERING } from '../../../modules/local/maf_rna_filtering/main'

workflow MAF_FILTERING_RNA {
    take:
    maf_to_filter             // maf from first pass
    maf_to_filter_realigned   // maf from realignment (second pass) [OPT?]
    fasta
    fasta_fai
    input_sample

    main:
    versions = Channel.empty()
    maf_to_filter.dump(tag:"maf_to_filter")
    maf_to_filter_realigned.dump(tag:"maf_to_filter_realigned")
    if (params.step in ['mapping', 'markduplicates', 'splitncigar',
    'prepare_recalibration', 'recalibrate', 'variant_calling', 'norm', 'consensus',
    'realignment', 'rna_filtering'] && (params.tools && params.tools.split(",").contains("rna_filtering"))) {

        if (params.step == 'rna_filtering') { maf_to_filter = input_sample} // TODO: not implemented yet
        else {
            if (params.step =="realignment"){
                maf_to_filter = maf_to_filter_realigned // there is no first maf
                maf_to_filter_realigned = null
            }
            if (maf_to_filter_realigned){
                // RNA filtering after realignment
                maf_to_cross_first_pass = maf_to_filter
                                            .map{meta, maf -> [meta.patient, meta, maf]}


                maf_to_cross_first_pass.dump(tag:"[STEP9: FILTERING] maf_to_cross_first_pass")

                maf_to_cross_realignment = maf_to_filter_realigned
                                            .map{meta, maf -> [meta.patient, meta, maf]}
                maf_to_cross_realignment.dump(tag:"[STEP9: FILTERING] maf_to_cross_realignment")
                maf_to_cross_first_pass
                                .cross(maf_to_cross_realignment).dump(tag:"[STEP9: FILTERING] maf_to_cross_crossed_pass")
                maf_crossed = maf_to_cross_first_pass
                                .cross(maf_to_cross_realignment)
                                .map{first, second ->
                                        def meta = [:]
                                        meta.patient    = first[0]
                                        meta.first_id   = first[1].id
                                        meta.second_id  = second[1].id
                                        meta.status     = first[1].status
                                        meta.tumor_id   = first[1].id.split("_vs_")[0]
                                        meta.id         = first[1].id.split("_vs_")[0] + "_with_" + second[1].id.split("_vs_")[0]
                                        meta.normal_id  = first[1].id.split("_vs_")[1]
                                        [meta, first[2], second[2]]
                                        }
            } else {
            maf_to_filter.dump(tag:"maf_to_filter2")
            maf_crossed = maf_to_filter.map{meta, maf -> [meta, maf, []]}

            }
        }
        maf_crossed.dump(tag:"[STEP9: FILTERING] maf_crossed")
//        maf_to_filter_status_dna = maf_to_filter_status.dna.map{meta, maf -> [meta, maf, file("NO_FILE.maf")]}
//        maf_to_filter_status_dna.dump(tag:"[STEP9: FILTERING] maf_to_filter_status_dna")
//        maf_crossed = maf_crossed.mix(maf_to_filter_status.dna)
        RNA_FILTERING(maf_crossed,
                    fasta,
                    fasta_fai)
        versions = versions.mix(RNA_FILTERING.out.versions)
    }


    emit:
        versions            = versions // channel: [ versions.yml ]



}
