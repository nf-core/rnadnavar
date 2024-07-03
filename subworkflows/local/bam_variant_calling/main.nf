//
// Variant Calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Variant calling on tumor/normal pair
include { BAM_VARIANT_CALLING_SOMATIC                 } from '../bam_variant_calling_somatic/main'
// QC on VCF files
include { VCF_QC_BCFTOOLS_VCFTOOLS                    } from '../vcf_qc_bcftools_vcftools/main'
// Create samplesheet to restart from different steps
include { CHANNEL_VARIANT_CALLING_CREATE_CSV          } from '../channel_variant_calling_create_csv/main'



workflow BAM_VARIANT_CALLING {

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
    intervals_bed_gz_tbi_combined
    pon
    pon_tbi
    input_sample
    realignment
    no_intervals

    main:
    reports   = Channel.empty()
    versions  = Channel.empty()
    if (tools || realignment) {
        if ((params.step == 'annotate' || params.step == 'norm') && (!realignment)) {
            cram_variant_calling_pair   = Channel.empty()
            vcf_to_normalise            = input_sample
            contamination_table_mutect2 = Channel.empty()
            segmentation_table_mutect2  = Channel.empty()
            artifact_priors_mutect2     = Channel.empty()

        } else {

            //
            // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
            //
            cram_variant_calling_status = cram_variant_calling.branch{
                normal: it[0].status == 0
                tumor:  it[0].status >= 1  // DNA and RNA should NOT have same sample id
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

            cram_variant_calling_pair.dump(tag:"cram_variant_calling_pair")

            // PAIR VARIANT CALLING
            BAM_VARIANT_CALLING_SOMATIC(
                tools,
                cram_variant_calling_pair,
                fasta,
                fasta_fai,
                dict,
                germline_resource,
                germline_resource_tbi,
                intervals,
                intervals_bed_gz_tbi,
                intervals_bed_combined,
                intervals_bed_gz_tbi_combined,
                pon,
                pon_tbi,
                params.joint_mutect2,
                realignment,
                no_intervals
                )


            // Gather vcf files for annotation and QC
            vcf_to_normalise            = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC.out.vcf_all)
            contamination_table_mutect2 = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC.out.contamination_table_mutect2)
            segmentation_table_mutect2  = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC.out.segmentation_table_mutect2)
            artifact_priors_mutect2     = Channel.empty().mix(BAM_VARIANT_CALLING_SOMATIC.out.artifact_priors_mutect2)

            // Gather used variant calling softwares versions
            versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC.out.versions)

            CHANNEL_VARIANT_CALLING_CREATE_CSV(vcf_to_normalise, "variantcalled")
        }

        // QC
        VCF_QC_BCFTOOLS_VCFTOOLS(vcf_to_normalise, intervals_bed_combined)
        versions = versions.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.versions)

        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.bcftools_stats.collect{ meta, stats -> stats })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_counts.collect{ meta, counts -> counts })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_qual.collect{ meta, qual -> qual })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_filter_summary.collect{ meta, summary -> summary })

    } else{

        cram_variant_calling_pair   = Channel.empty()
        vcf_to_normalise            = Channel.empty()
        contamination_table_mutect2 = Channel.empty()
        segmentation_table_mutect2  = Channel.empty()
        artifact_priors_mutect2     = Channel.empty()

    }

    emit:
    cram_variant_calling_pair        = cram_variant_calling_pair
    vcf_to_normalise                 = vcf_to_normalise
    contamination_table              = contamination_table_mutect2
    segmentation_table               = segmentation_table_mutect2
    artifact_priors                  = artifact_priors_mutect2
    reports                          = reports
    versions                         = versions

}
