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
    'prepare_recalibration', 'recalibrate', 'variant_calling', 'norm', 'consensus', 'filtering',
    'realignment', 'rna_filtering'] && (params.tools && params.tools.split(",").contains("rna_filtering"))) {

        // Allow direct restart from MAF input when RNA-specific filtering is the entry step.
        // If both first-pass and realigned MAFs are provided, distinguish them here using
        // the historical `_realign` suffix so the same crossing logic can be reused.
        if (params.step == 'rna_filtering') {
            input_sample_type = input_sample.branch { meta, maf ->
                realigned: meta.id.contains('_realign')
                first_pass: true
            }
            maf_to_filter = input_sample_type.first_pass
            maf_to_filter_realigned = input_sample_type.realigned
        } else {
            if (params.step == "realignment"){
                maf_to_filter = maf_to_filter_realigned // there is no first maf
                maf_to_filter_realigned = Channel.empty()
            }
        }

        // Pair first-pass and realigned MAFs locally by patient. Using a keyed join with
        // remainder avoids probing channel truthiness and still supports both restart modes:
        // with a realigned MAF we emit [meta, first_maf, realigned_maf], otherwise
        // [meta, first_maf, []] for single-input RNA filtering.
        maf_to_cross_first_pass = maf_to_filter
                                    .map{meta, maf -> [meta.patient, meta, maf]}
        maf_to_cross_realignment = maf_to_filter_realigned
                                    .map{meta, maf -> [meta.patient, meta, maf]}

        maf_to_cross_first_pass.dump(tag:"maf_to_cross_first_pass")
        maf_to_cross_realignment.dump(tag:"maf_to_cross_realignment")

        maf_crossed = maf_to_cross_first_pass
                        .join(maf_to_cross_realignment, remainder: true)
                        .map{ patient, first_meta, first_maf, second_meta, second_maf ->
                                if (second_meta) {
                                    def meta = [:]
                                    meta.patient    = patient
                                    meta.first_id   = first_meta.id
                                    meta.second_id  = second_meta.id
                                    meta.status     = first_meta.status
                                    meta.tumor_id   = first_meta.id.split("_vs_")[0]
                                    meta.id         = first_meta.id.split("_vs_")[0] + "_with_" + second_meta.id.split("_vs_")[0]
                                    meta.normal_id  = first_meta.id.split("_vs_")[1]
                                    [meta, first_maf, second_maf]
                                } else {
                                    [first_meta, first_maf, []]
                                }
                            }
        maf_crossed.dump(tag:"[STEP9: FILTERING] maf_crossed")
//        maf_to_filter_status_dna = maf_to_filter_status.dna.map{meta, maf -> [meta, maf, file("NO_FILE.maf")]}
//        maf_to_filter_status_dna.dump(tag:"[STEP9: FILTERING] maf_to_filter_status_dna")
//        maf_crossed = maf_crossed.mix(maf_to_filter_status.dna)
        RNA_FILTERING(maf_crossed,
                    fasta)
        versions = versions.mix(RNA_FILTERING.out.versions)
    }


    emit:
        versions            = versions // channel: [ versions.yml ]



}
