include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM_VARIANTCALLING } from '../../modules/nf-core/modules/samtools/convert/main'
include { PAIR_VARIANT_CALLING                                  } from './pair_variant_calling'
include { VCF_QC                                               } from '../nf-core/vcf_qc'
include { VARIANTCALLING_CSV                                   } from './variantcalling_csv'


workflow VARIANT_CALLING {

    take:
        tools
        ch_cram_variant_calling
        fasta
        fasta_fai
        dbsnp
        dbsnp_tbi
        dict
        germline_resource
        germline_resource_tbi
        intervals
        intervals_bed_gz_tbi
        intervals_bed_combined
        pon
        pon_tbi
        ch_input_sample

    main:
        ch_reports   = Channel.empty()
        ch_versions  = Channel.empty()

        // get input for variant calling
        if (params.step == 'variant_calling') {
            // if input is a BAM file for variant calling we need to convert it to CRAM
            ch_input_sample.branch{
                    bam: it[0].data_type == "bam"
                    cram: it[0].data_type == "cram"
                }.set{ch_convert}
            //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            SAMTOOLS_BAMTOCRAM_VARIANTCALLING(ch_convert.bam, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM_VARIANTCALLING.out.versions)
            ch_cram_variant_calling = Channel.empty().mix(SAMTOOLS_BAMTOCRAM_VARIANTCALLING.out.alignment_index, ch_convert.cram)
        }

        if (params.step == 'annotate') {
            // no variant calling will be performed
            ch_cram_variant_calling = Channel.empty()
        }
        // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
        ch_cram_variant_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status >= 1
        }.set{ch_cram_variant_calling_status}

        // All Germline samples -- will be the same for DNA and RNA
        ch_cram_variant_calling_normal_to_cross = ch_cram_variant_calling_status.normal.map{ meta, cram, crai -> [meta.patient, meta, cram, crai] }
//            ch_cram_variant_calling_normal_to_cross.dump(tag: "[STEP3_VARIANTCALLING] normals to cross")
        // All tumor samples
        ch_cram_variant_calling_tumor_pair_to_cross = ch_cram_variant_calling_status.tumor.map{ meta, cram, crai -> [meta.patient, meta, cram, crai] }
//            ch_cram_variant_calling_tumor_pair_to_cross.dump(tag: "[STEP3_VARIANTCALLING] tumors to cross")

        // Tumor - normal pairs
        // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
        ch_cram_variant_calling_pair = ch_cram_variant_calling_normal_to_cross.cross(ch_cram_variant_calling_tumor_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]
                meta.patient    = normal[0]
                meta.normal_id  = normal[1].sample
                meta.tumor_id   = tumor[1].sample
                meta.status     = tumor[1].status
                meta.id         = "${meta.tumor_id}_vs_${meta.normal_id}".toString()
                meta.alleles    = null
                [meta, normal[2], normal[3], tumor[2], tumor[3]]
            }
//            ch_cram_variant_calling_pair.dump(tag:"[STEP3_VARIANTCALLING] variant_calling_pairs")
        PAIR_VARIANT_CALLING(
            params.tools,
            ch_cram_variant_calling_pair,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            germline_resource,
            germline_resource_tbi,
            intervals,
            intervals_bed_gz_tbi,
            intervals_bed_combined,
            pon,
            pon_tbi,
            params.highconfidence,
            params.actionablepanel,
            params.knownhot,
            params.ensbl_sage,
            params.skip_tools
        )
        ch_versions = ch_versions.mix(PAIR_VARIANT_CALLING.out.versions)

        // Gather vcf files for annotation and QC
        vcf_to_normalize = Channel.empty()
        vcf_to_normalize = vcf_to_normalize.mix(PAIR_VARIANT_CALLING.out.strelka_vcf)
        vcf_to_normalize = vcf_to_normalize.mix(PAIR_VARIANT_CALLING.out.mutect2_vcf)
        vcf_to_normalize = vcf_to_normalize.mix(PAIR_VARIANT_CALLING.out.freebayes_vcf)
        vcf_to_normalize = vcf_to_normalize.mix(PAIR_VARIANT_CALLING.out.sage_vcf)
    //        ch_cram_variant_calling_pair.dump(tag:"[STEP3_VARIANTCALLING] all_vcfs")

        //QC
        if (tools.split(',').contains('vcf_qc')) {
            VCF_QC(vcf_to_normalize, intervals_bed_combined)
            VARIANTCALLING_CSV(vcf_to_normalize)

            ch_versions = ch_versions.mix(VCF_QC.out.versions)
            ch_reports  = ch_reports.mix(VCF_QC.out.bcftools_stats.collect{meta, stats -> stats})
            ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_counts.collect{ meta, counts -> counts})
            ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_tstv_qual.collect{ meta, qual -> qual })
            ch_reports  = ch_reports.mix(VCF_QC.out.vcftools_filter_summary.collect{meta, summary -> summary})
        }


    emit:
        cram_vc_pair        = ch_cram_variant_calling_pair
        vcf                 = vcf_to_normalize
        contamination_table = PAIR_VARIANT_CALLING.out.contamination_table
        segmentation_table  = PAIR_VARIANT_CALLING.out.segmentation_table
        artifact_priors     = PAIR_VARIANT_CALLING.out.artifact_priors
        reports             = ch_reports
        versions            = ch_versions

}