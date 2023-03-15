//
// NORMALIZATION OF VCF VARIANTS WITH VT
//

include { VT_DECOMPOSE } from '../../modules/local/vt/decompose/main'
include { VT_NORMALIZE } from '../../modules/local/vt/normalize/main'

workflow NORMALIZE {
    take:
        vcf
        fasta
        ch_input_sample

    main:
        ch_versions = Channel.empty()

        if (params.step in ['mapping', 'markduplicates', 'splitncigar', 'prepare_recalibration', 'recalibrate', 'variant_calling', 'normalize'] ) {
            if (params.step == 'normalize')
            // TODO - test this
                vcf_to_normalize = ch_input_sample
            else {
                vcf_to_normalize = vcf
                }
            vcf_to_normalize.dump(tag:'[STEP4] vcf_to_normalize')
            ch_vcf_decomp  = Channel.empty()
            ch_vcf_norm  = Channel.empty()
            // Separate variants
            VT_DECOMPOSE(vcf)

            ch_vcf_decomp = ch_vcf_decomp.mix(VT_DECOMPOSE.out.vcf)
            ch_versions = ch_versions.mix(VT_DECOMPOSE.out.versions.first())

            // Normalize variants
            VT_NORMALIZE(ch_vcf_decomp,
                         fasta)

            ch_vcf_norm = ch_vcf_norm.mix(VT_NORMALIZE.out.vcf)
            ch_versions = ch_versions.mix(VT_NORMALIZE.out.versions.first())
            ch_vcf_norm.dump(tag:'[STEP4] vcf_normalized')
        }

    emit:
        vcf         = ch_vcf_norm // channel: [ [meta], vcf ]
        versions    = ch_versions // channel: [ versions.yml ]

}
