

include { DNA_FILTERING                                        } from '../../../modules/local/dna_filtering'



workflow FILTER_DNA {

    take:
        maf
        whitelist

    main:
        ch_reports   = Channel.empty()
        ch_versions  = Channel.empty()
        DNA_FILTERING(maf, whitelist)

    emit:
        filter_maf = DNA_FILTERING.out.maf
}