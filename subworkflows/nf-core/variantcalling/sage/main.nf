include { BCFTOOLS_SORT                                } from '../../../../modules/nf-core/modules/bcftools/sort/main'
include { GATK4_MERGEVCFS as MERGE_SAGE           } from '../../../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { SAGE                                    } from '../../../../modules/nf-core/modules/sage/main'
include { TABIX_TABIX as TABIX_VC_SAGE            } from '../../../../modules/nf-core/modules/tabix/tabix/main'
//include { TABIX_BGZIPTABIX as BGZIPTABIX_VC_SAGE            } from '../../../../modules/nf-core/modules/tabix/bgziptabix/main'


workflow RUN_SAGE {
    take:
        cram                     // channel: [mandatory] [meta, cram, crai, [], [], interval]
        dict
        fasta                    // channel: [mandatory]
        fasta_fai
        highconfidence
        actionablepanel
        knownhot
        ensbl_sage

    main:

        ch_versions = Channel.empty()


        SAGE(
            cram,
            fasta,
            fasta_fai,
            dict,
            highconfidence,
            actionablepanel,
            knownhot,
            ensbl_sage)

        BCFTOOLS_SORT(SAGE.out.vcf)
        BCFTOOLS_SORT.out.vcf.branch{
                intervals:    it[0].num_intervals > 1
                no_intervals: it[0].num_intervals <= 1
            }.set{bcftools_vcf_out}

        // Only when no intervals
        TABIX_VC_SAGE(bcftools_vcf_out.no_intervals)

        // Only when using intervals
        MERGE_SAGE(
            bcftools_vcf_out.intervals
                .map{ meta, vcf ->

                    new_meta = meta.tumor_id ? [
                                                    id:             meta.tumor_id + "_vs_" + meta.normal_id,
                                                    normal_id:      meta.normal_id,
                                                    num_intervals:  meta.num_intervals,
                                                    patient:        meta.patient,
                                                    status:         meta.status,
                                                    tumor_id:       meta.tumor_id,
                                                    alleles:        meta.alleles
                                                ]
                                            :   [
                                                    id:             meta.sample,
                                                    num_intervals:  meta.num_intervals,
                                                    patient:        meta.patient,
                                                    sample:         meta.sample,
                                                    status:         meta.status,
                                                    alleles:        meta.alleles
                                                ]
                    [groupKey(new_meta, meta.num_intervals), vcf]
                }.groupTuple(),
            dict
        )

        // Mix output channels for "no intervals" and "with intervals" results
        sage_vcf = Channel.empty().mix(
                            MERGE_SAGE.out.vcf,
                            bcftools_vcf_out.no_intervals)
                        .map{ meta, vcf ->
                            [ [
                                id:             meta.id,
                                normal_id:      meta.normal_id,
                                num_intervals:  meta.num_intervals,
                                patient:        meta.patient,
                                status:         meta.status,
                                tumor_id:       meta.tumor_id,
                                alleles:        meta.alleles,
                                variantcaller:  "sage"
                            ],
                                vcf]
                        }


        ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
        ch_versions = ch_versions.mix(MERGE_SAGE.out.versions)
        ch_versions = ch_versions.mix(SAGE.out.versions)

    emit:
        sage_vcf = sage_vcf
        versions = ch_versions
}
