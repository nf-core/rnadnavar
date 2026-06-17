//
// Normalise VCFs with VT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// VT steps
include { VT_DECOMPOSE                        } from '../../../modules/nf-core/vt/decompose/main'
include { VT_NORMALIZE                        } from '../../../modules/nf-core/vt/normalize/main'
include { HTSLIB_BGZIPTABIX as INDEX_VCF_DECOMPOSED } from '../../../modules/nf-core/htslib/bgziptabix/main'
// Create samplesheet to restart from different steps
include { CHANNEL_VARIANT_CALLING_CREATE_CSV  } from '../channel_variant_calling_create_csv/main'


workflow VCF_NORMALIZE {
    take:
    vcf_to_normalize
    fasta
    fasta_fai
    input_sample
    realignment

    main:
    version          = Channel.empty()

    if (params.step == 'norm') vcf_to_normalize = input_sample

    if ((params.step in ['mapping', 'markduplicates', 'splitncigar',
                        'prepare_recalibration', 'recalibrate',
                        'variant_calling', 'norm'] &&
                        ((params.tools && params.tools.split(",").contains("consensus")))) ||
                        realignment) {

        vcf_decomposed  = Channel.empty()
        vcf_to_normalize = vcf_to_normalize.map{meta, vcf -> [meta, vcf, []]} // vt accepts intervals, not in use for now
        // Separate variants
        VT_DECOMPOSE(vcf_to_normalize)

        vcf_decomposed = vcf_decomposed.mix(VT_DECOMPOSE.out.vcf)
        version = version.mix(VT_DECOMPOSE.out.versions_vt)

        // VT normalize now expects a compressed VCF together with its tabix index.
        INDEX_VCF_DECOMPOSED(
            vcf_decomposed.map { meta, vcf_file -> [meta, vcf_file, [], []] },
            'compress',
            true,
            'vcf'
        )
        version = version.mix(INDEX_VCF_DECOMPOSED.out.versions_htslib)
        version = version.mix(INDEX_VCF_DECOMPOSED.out.versions_xz)

        // Normalise variants
        vcf_decomposed_indexed = INDEX_VCF_DECOMPOSED.out.output
            .join(INDEX_VCF_DECOMPOSED.out.index, failOnMismatch: true)
            .map { meta, vcf_file, tbi_file -> [meta, vcf_file, tbi_file, []] }

        VT_NORMALIZE(vcf_decomposed_indexed,
                    fasta, fasta_fai) // fai not necessary?

        vcf_to_consensus = VT_NORMALIZE.out.vcf
        version = version.mix(VT_NORMALIZE.out.versions_vt)

        CHANNEL_VARIANT_CALLING_CREATE_CSV(vcf_to_consensus, "normalized")

    } else {
        vcf_to_consensus = vcf_to_normalize
    }

    emit:
    vcf         = vcf_to_consensus // channel: [ [meta], vcf ]
    versions    = version // channel: [ versions.yml ]

}
