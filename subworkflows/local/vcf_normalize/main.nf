//
// Normalise VCFs with VT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// VT steps
include { VT_DECOMPOSE                        } from '../../../modules/nf-core/vt/decompose/main'
include { VT_NORMALIZE                        } from '../../../modules/nf-core/vt/normalize/main'
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

        // Normalise variants
        vcf_decomposed = vcf_decomposed.map{meta,vcf -> [meta, vcf, [],[]]} // tbi not necessary, vt accepts intervals, not in use for now
        VT_NORMALIZE(vcf_decomposed,
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
