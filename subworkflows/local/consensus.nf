//
// Make a consensus out of different VCF files
//

include { RUN_CONSENSUS } from '../../modules/local/run_consensus'


workflow CONSENSUS {
    take:
        vcfs
        callers

    main:
        ch_versions = Channel.empty()
        consensus_vcf = Channel.empty()
        ch_versions = Channel.empty()
        vcfs.dump(tag:'vcfs in consensus')
        callers.dump(tag:'callers in consensus')
        RUN_CONSENSUS ( vcfs, callers )

     //   consensus_vcf = RUN_CONSENSUS.out.vcf

     //   ch_versions = ch_versions.mix(RUN_CONSENSUS.out.versions.first())

    emit:
        vcf    = consensus_vcf // channel: [ [meta], vcf ]
        versions    = ch_versions // channel: [ versions.yml ]
}