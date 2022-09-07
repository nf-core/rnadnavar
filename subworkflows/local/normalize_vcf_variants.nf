//
// NORMALIZATION OF VCF VARIANTS WITH VT
//

include { VT_DECOMPOSE } from '../../modules/local/vt/decompose/main'
include { VT_NORMALIZE } from '../../modules/local/vt/normalize/main'

workflow NORMALIZE_VCF {
    take:
    vcf
    fasta

    main:

    ch_versions = Channel.empty()

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

    emit:
        vcf         = ch_vcf_norm // channel: [ [meta], vcf ]
        versions    = ch_versions // channel: [ versions.yml ]

}
