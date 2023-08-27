include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM_VARIANTCALLING } from '../../modules/nf-core/modules/samtools/convert/main'
include { PAIR_VARIANT_CALLING                                  } from './pair_variant_calling'
include { VCF_QC                                               } from '../nf-core/vcf_qc'
include { VARIANTCALLING_CSV                                   } from './variantcalling_csv'


workflow VARIANT_CALLING {

    take:
    tools
    cram_variant_calling
    fasta
    fasta_fai
    dict
    germline_resource
    germline_resource_tbi
    intervals
    intervals_bed_gz_tbi
    intervals_bed_combined
    pon
    pon_tbi
    input_sample

    main:
    reports   = Channel.empty()
    versions  = Channel.empty()

    if (params.step == 'variant_calling') {

        input_variant_calling_convert = input_sample.branch{
            bam:  it[0].data_type == "bam"
            cram: it[0].data_type == "cram"
        }

        // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
        BAM_TO_CRAM(input_variant_calling_convert.bam, fasta, fasta_fai)
        versions = versions.mix(BAM_TO_CRAM.out.versions)

        cram_variant_calling = Channel.empty().mix(BAM_TO_CRAM.out.alignment_index, input_variant_calling_convert.cram)

    }

    if (params.tools) {
        if (params.step == 'annotate') cram_variant_calling = Channel.empty()

	    //
        // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
        //
        cram_variant_calling_status = cram_variant_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
            rna:    it[0].status == 2
        }

            // All Germline samples
        cram_variant_calling_normal_to_cross = cram_variant_calling_status.normal.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // All tumor samples
        cram_variant_calling_pair_to_cross = cram_variant_calling_status.tumor.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // Tumor only samples
        // 1. Group together all tumor samples by patient ID [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ]

        // Downside: this only works by waiting for all tumor samples to finish preprocessing, since no group size is provided
        cram_variant_calling_tumor_grouped = cram_variant_calling_pair_to_cross.groupTuple()

        // 2. Join with normal samples, in each channel there is one key per patient now. Patients without matched normal end up with: [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ], null ]
        cram_variant_calling_tumor_joined = cram_variant_calling_tumor_grouped.join(cram_variant_calling_normal_to_cross, failOnDuplicate: true, remainder: true)

        // 3. Filter out entries with last entry null
        cram_variant_calling_tumor_filtered = cram_variant_calling_tumor_joined.filter{ it ->  !(it.last()) }

        // 4. Transpose [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ] back to [ patient1, meta1, [ cram1, crai1 ], null ] [ patient1, meta2, [ cram2, crai2 ], null ]
        // and remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ]
        cram_variant_calling_tumor_only = cram_variant_calling_tumor_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

		only_paired_variant_calling = true // for now only this supported
	    if (only_paired_variant_calling) {
            // Normal only samples

            // 1. Join with tumor samples, in each channel there is one key per patient now. Patients without matched tumor end up with: [ patient1, [ meta1 ], [ cram1, crai1 ], null ] as there is only one matched normal possible
            cram_variant_calling_normal_joined = cram_variant_calling_normal_to_cross.join(cram_variant_calling_tumor_grouped, failOnDuplicate: true, remainder: true)

            // 2. Filter out entries with last entry null
            cram_variant_calling_normal_filtered = cram_variant_calling_normal_joined.filter{ it ->  !(it.last()) }

            // 3. Remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ] (no transposing needed since only one normal per patient ID)
            cram_variant_calling_status_normal = cram_variant_calling_normal_filtered.map{ it -> [it[1], it[2], it[3]] }
		} else {
            cram_variant_calling_status_normal = cram_variant_calling_status.normal
        }

        // Tumor - normal pairs
        // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
        cram_variant_calling_pair = cram_variant_calling_normal_to_cross.cross(cram_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.status     = tumor[1].status
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
            }


	    // PAIR VARIANT CALLING
	    BAM_VARIANT_CALLING_SOMATIC(
	        tools,
	        cram_variant_calling_pair,
	        dict,
	        fasta,
	        fasta_fai,
	        germline_resource,
	        germline_resource_tbi,
	        intervals,
	        intervals_bed_gz_tbi,
	        intervals_bed_combined,
	        pon,
	        pon_tbi
	        )

	    // POST VARIANTCALLING
        POST_VARIANTCALLING(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_all,
                            params.concatenate_vcfs)

	    // Gather vcf files for annotation and QC
	    vcf_to_normalize = Channel.empty()
	    vcf_to_normalize = vcf_to_normalize.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.vcf_all)

	    // QC
	    VCF_QC_BCFTOOLS_VCFTOOLS(vcf_to_normalize, intervals_bed_combined)

	    reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.bcftools_stats.collect{ meta, stats -> stats })
	    reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_counts.collect{ meta, counts -> counts })
	    reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_qual.collect{ meta, qual -> qual })
	    reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_filter_summary.collect{ meta, summary -> summary })

	    CHANNEL_VARIANT_CALLING_CREATE_CSV(vcf_to_annotate)

		// Gather used variant calling softwares versions
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.versions)
        versions = versions.mix(POST_VARIANTCALLING.out.versions)
        versions = versions.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.versions)
    }


    emit:
    cram_vc_pair        = ch_cram_variant_calling_pair
    vcf                 = vcf_to_normalize
    contamination_table = PAIR_VARIANT_CALLING.out.contamination_table
    segmentation_table  = PAIR_VARIANT_CALLING.out.segmentation_table
    artifact_priors     = PAIR_VARIANT_CALLING.out.artifact_priors
    reports             = reports
    versions            = versions

}