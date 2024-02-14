//
// GATK pre-processing best practices
//
// Markduplicates
include { SAMTOOLS_CONVERT as BAM_TO_CRAM                      } from '../../../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM                      } from '../../../modules/nf-core/samtools/convert/main'
include { BAM_MARKDUPLICATES                                   } from '../../local/bam_markduplicates/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV                    } from '../../local/channel_markduplicates_create_csv/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD           } from '../../local/cram_qc_mosdepth_samtools/main'
// Splitncigarreads
include { BAM_SPLITNCIGARREADS                                 } from '../../local/bam_splitncigarreads/main'
include { CHANNEL_SPLITNCIGARREADS_CREATE_CSV                  } from '../../local/channel_splitncigarreads_create_csv/main'
// Create recalibration tables
include { BAM_BASERECALIBRATOR                                 } from '../../local/bam_baserecalibrator/main'
include { CHANNEL_BASERECALIBRATOR_CREATE_CSV                  } from '../../local/channel_baserecalibrator_create_csv/main'
// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR                                        } from '../../local/bam_applybqsr/main'
include { CHANNEL_APPLYBQSR_CREATE_CSV                         } from '../../local/channel_applybqsr_create_csv/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL           } from '../../local/cram_qc_mosdepth_samtools/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL                } from '../../../modules/nf-core/samtools/convert/main'


workflow BAM_GATK_PREPROCESSING {
    take:
    input_sample                  // channel: [optional]  input from CSV if applicable
    bam_mapped                    // channel: [mandatory] bam_mapped
    cram_mapped                   // channel: [mandatory] cram_mapped
    fasta                         // channel: [mandatory] fasta
    fasta_fai                     // channel: [mandatory] fasta_fai
    dict                          // channel: [mandatory] dict
    known_sites_indels            // channel: [optional]  known_sites
    known_sites_indels_tbi        // channel: [optional]  known_sites
    germline_resource             // channel: [optional]  germline_resource
    germline_resource_tbi         // channel: [optional]  germline_resource_tbi
    intervals                     // channel: [mandatory] intervals/target regions
    intervals_for_preprocessing   // channel: [mandatory] intervals/wes
    intervals_and_num_intervals   // channel: [mandatory] [ intervals, num_intervals ] (or [ [], 0 ] if no intervals)
    second_run                    // boolean


    main:
    reports   = Channel.empty()
    versions  = Channel.empty()
    cram_variant_calling = Channel.empty()
    // Markduplicates
    if (params.step in ['mapping', 'markduplicates'] || second_run) {

        cram_markduplicates_no_spark = Channel.empty()
        // make sure data_type is present
        bam_mapped = bam_mapped.map { infoMap ->
                                    infoMap[0] = (infoMap[0] + [data_type: "bam"])
                                    infoMap }
        cram_mapped = cram_mapped.map { infoMap ->
                                    infoMap[0] = (infoMap[0] + [data_type: "cram"])
                                    infoMap }
        // cram_for_markduplicates will contain bam mapped with FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        cram_for_markduplicates = params.step == 'mapping' || second_run ? bam_mapped : input_sample.map{ meta, input, index -> [ meta, input ] }

        // if no MD is done, then run QC on mapped & converted CRAM files
        // or the input BAM (+converted) or CRAM files
        cram_skip_markduplicates = Channel.empty()

        // Should it be possible to restart from converted crams?
        // For now, conversion from bam to cram is only done when skipping markduplicates
        if ((params.skip_tools && params.skip_tools.split(',').contains('markduplicates')) && !second_run) {
            if (params.step == 'mapping') {
                cram_skip_markduplicates = cram_for_markduplicates
            } else {
                input_markduplicates_convert = cram_for_markduplicates.branch{
                    bam:  it[0].data_type == "bam"
                    cram: it[0].data_type == "cram"
                }

                // Convert any input BAMs to CRAM
                BAM_TO_CRAM(input_markduplicates_convert.bam, fasta, fasta_fai)
                versions = versions.mix(BAM_TO_CRAM.out.versions)

                cram_skip_markduplicates = Channel.empty().mix(input_markduplicates_convert.cram, BAM_TO_CRAM.out.alignment_index)
            }

            CRAM_QC_NO_MD(cram_skip_markduplicates, fasta, intervals_for_preprocessing)

            // Gather QC reports
            reports = reports.mix(CRAM_QC_NO_MD.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(CRAM_QC_NO_MD.out.versions)
        } else {
            cram_for_markduplicates.dump(tag:"cram_for_markduplicates")
            BAM_MARKDUPLICATES(
                cram_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(BAM_MARKDUPLICATES.out.versions)
        }

        // ch_md_cram_for_restart contains crams from markduplicates
        ch_md_cram_for_restart = Channel.empty().mix(cram_markduplicates_no_spark)
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }
        // If params.save_output_as_bam, then convert CRAM files to BAM
        CRAM_TO_BAM(ch_md_cram_for_restart, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM.out.versions)

        // CSV should be written for the file actually out, either CRAM or BAM
        // Create CSV to restart from this step
        csv_subfolder = 'markduplicates'
        params.save_output_as_bam ? CHANNEL_MARKDUPLICATES_CREATE_CSV(CRAM_TO_BAM.out.alignment_index, csv_subfolder, params.outdir, params.save_output_as_bam) : CHANNEL_MARKDUPLICATES_CREATE_CSV(ch_md_cram_for_restart, csv_subfolder, params.outdir, params.save_output_as_bam)
    } else {
        ch_md_cram_for_restart   = Channel.empty().mix(input_sample)
        cram_skip_markduplicates = Channel.empty().mix(input_sample)
    }


    // SplitNCigarReads for RNA
    cram_skip_splitncigar = ch_md_cram_for_restart?: cram_skip_markduplicates
    if (params.step in ['mapping', 'markduplicates', 'splitncigar'] || second_run) {
        if (params.step == "splitncigarreads"){
                // Support if starting from BAM or CRAM files
                input_sncr_convert = input_sample.branch{
                    bam:  it[0].data_type == "bam"
                    cram: it[0].data_type == "cram"
                }

                input_sncr_convert  = input_sncr_convert.bam.map{ meta, bam, bai, table -> [ meta, bam, bai ] }
                // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
                BAM_TO_CRAM(input_sncr_convert, fasta, fasta_fai)
                versions = versions.mix(BAM_TO_CRAM.out.versions)

                sncr_cram_from_bam = BAM_TO_CRAM.out.alignment_index
                    // Make sure correct data types are carried through
                    .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

                cram_skip_splitncigar = Channel.empty().mix(sncr_cram_from_bam, input_sncr_convert.cram)
        }

        // cram_for_bam_splitncigar contains either:
        // - crams from markduplicates
        // - crams converted from bam mapped when skipping markduplicates
        // - input cram files, when start from step markduplicates
        cram_skip_splitncigar = cram_skip_splitncigar
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }
        cram_splitncigar_no_spark = Channel.empty()
        if (!(params.skip_tools && params.skip_tools.split(',').contains('splitncigar')) || second_run) {
            cram_skip_splitncigar.dump(tag:"cram_skip_splitncigar")
            cram_for_splitncigar_status = cram_skip_splitncigar.branch{
                                                dna:  it[0].status < 2
                                                rna:  it[0].status >= 2
                                            }
            BAM_SPLITNCIGARREADS (
                cram_for_splitncigar_status.rna,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals
            )

            cram_splitncigar_no_spark = BAM_SPLITNCIGARREADS.out.cram.mix(cram_for_splitncigar_status.dna)

            // Gather used softwares versions
            versions = versions.mix(BAM_SPLITNCIGARREADS.out.versions)

            cram_skip_splitncigar = Channel.empty()
        }

        // cram_splitncigar_no_spark contains crams from splitncigar or markduplicates if applicable
        ch_sncr_cram_for_restart = Channel.empty().mix(cram_splitncigar_no_spark).mix(cram_skip_splitncigar)
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }
    } else {
        ch_sncr_cram_for_restart = Channel.empty()
    }

    // BQSR
    ch_cram_for_bam_baserecalibrator = ch_sncr_cram_for_restart?: cram_skip_splitncigar
    if (params.step in ['mapping', 'markduplicates', 'splitncigar', 'prepare_recalibration'] & !second_run) {

        // Run if starting from step "prepare_recalibration". This will not run for second pass
        if (params.step == 'prepare_recalibration' && !second_run) {

            // Support if starting from BAM or CRAM files
            input_prepare_recal_convert = input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }

            // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            BAM_TO_CRAM(input_prepare_recal_convert.bam, fasta, fasta_fai)
            versions = versions.mix(BAM_TO_CRAM.out.versions)

            sncr_cram_from_bam = BAM_TO_CRAM.out.alignment_index
                // Make sure correct data types are carried through
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_cram_for_bam_baserecalibrator, input_prepare_recal_convert.cram)
            ch_sncr_cram_for_restart = sncr_cram_from_bam

        } else {
            // ch_cram_for_bam_baserecalibrator contains either:
            // - crams from markduplicates
            // - crams from splitncigarreads
            // - crams converted from bam mapped when skipping markduplicates
            // - crams converted from bam mapped when skipping splitncigarreads
            // - input cram files, when start from step markduplicates
            // - input cram files, when start from step splitncigarreads
            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_sncr_cram_for_restart, cram_skip_splitncigar )
                // Make sure correct data types are carried through
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }
        }

        // Create recalibration tables
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator')) && !second_run) {

            ch_table_bqsr_no_spark = Channel.empty()
            ch_cram_for_bam_baserecalibrator.dump(tag:"ch_cram_for_bam_baserecalibrator")
            BAM_BASERECALIBRATOR(
                    ch_cram_for_bam_baserecalibrator,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals,
                    known_sites_indels,
                    known_sites_indels_tbi)

            ch_table_bqsr_no_spark = BAM_BASERECALIBRATOR.out.table_bqsr

            // Gather used softwares versions
            versions = versions.mix(BAM_BASERECALIBRATOR.out.versions)


            // ch_table_bqsr contains either:
            // - bqsr table from baserecalibrator
            ch_table_bqsr = Channel.empty().mix(
                ch_table_bqsr_no_spark)

            reports = reports.mix(ch_table_bqsr.collect{ meta, table -> table })
            cram_applybqsr = ch_cram_for_bam_baserecalibrator.join(ch_table_bqsr, failOnDuplicate: true, failOnMismatch: true)
            cram_applybqsr = cram_applybqsr.dump(tag:"cram_applybqsr")
            // Create CSV to restart from this step
            CHANNEL_BASERECALIBRATOR_CREATE_CSV(ch_sncr_cram_for_restart.join(ch_table_bqsr, failOnDuplicate: true), params.tools, params.skip_tools, params.save_output_as_bam, params.outdir)
        }
    } else{
        cram_variant_calling = ch_sncr_cram_for_restart
        cram_applybqsr       = ch_sncr_cram_for_restart

    }

    if (params.step in ['mapping', 'markduplicates', 'splitncigar','prepare_recalibration', 'recalibrate'] && !second_run) {

        // Run if starting from step "prepare_recalibration"
        if (params.step == 'recalibrate' && !second_run) {

            // Support if starting from BAM or CRAM files
            input_recal_convert = input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }

            // If BAM file, split up table and mapped file to convert BAM to CRAM
            input_only_table = input_recal_convert.bam.map{ meta, bam, bai, table -> [ meta, table ] }
            input_only_bam   = input_recal_convert.bam.map{ meta, bam, bai, table -> [ meta, bam, bai ] }

            // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            BAM_TO_CRAM(input_only_bam, fasta, fasta_fai)
            versions = versions.mix(BAM_TO_CRAM.out.versions)

            cram_applybqsr = Channel.empty().mix(
                BAM_TO_CRAM.out.alignment_index.join(input_only_table, failOnDuplicate: true, failOnMismatch: true),
                input_recal_convert.cram)
                // Join together converted cram with input tables
                .map{ meta, cram, crai, table -> [ meta + [data_type: "cram"], cram, crai, table ]}
        }
        if (!(params.skip_tools && params.skip_tools.split(',').contains('baserecalibrator')) && !second_run) {

            cram_variant_calling_no_spark = Channel.empty()
            cram_applybqsr.dump(tag:"cram_applybqsr")
            BAM_APPLYBQSR(
                cram_applybqsr,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals)

            cram_variant_calling_no_spark = BAM_APPLYBQSR.out.cram

            // Gather used softwares versions
            versions = versions.mix(BAM_APPLYBQSR.out.versions)

            cram_variant_calling = Channel.empty().mix(
                cram_variant_calling_no_spark)

            CRAM_QC_RECAL(
                cram_variant_calling,
                fasta,
                intervals_for_preprocessing)

            // Gather QC reports
            reports = reports.mix(CRAM_QC_RECAL.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(CRAM_QC_RECAL.out.versions)

            // If params.save_output_as_bam, then convert CRAM files to BAM
            CRAM_TO_BAM_RECAL(cram_variant_calling, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM_RECAL.out.versions)

            // CSV should be written for the file actually out out, either CRAM or BAM
            csv_recalibration = Channel.empty()
            csv_recalibration = params.save_output_as_bam ?  CRAM_TO_BAM_RECAL.out.alignment_index : cram_variant_calling

            // Create CSV to restart from this step
            CHANNEL_APPLYBQSR_CREATE_CSV(csv_recalibration)

        } else if (params.step == 'recalibrate') {
            // cram_variant_calling contains either:
            // - input bams converted to crams, if started from step recal + skip BQSR
            // - input crams if started from step recal + skip BQSR
            cram_variant_calling = Channel.empty().mix(
                BAM_TO_CRAM.out.alignment_index,
                input_recal_convert.cram.map{ meta, cram, crai, table -> [ meta, cram, crai ] })
        } else {
            // cram_variant_calling contains either:
            // - crams from markduplicates = ch_cram_for_bam_baserecalibrator if skip BQSR but not started from step recalibration
            cram_variant_calling = Channel.empty().mix(ch_cram_for_bam_baserecalibrator)
        }
    } else{
        cram_variant_calling = cram_applybqsr
    }


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
    // Remove lane from id (which is sample)
    cram_variant_calling = cram_variant_calling.map{meta, cram, crai ->
//													meta['read_group'] = meta['read_group'].replaceFirst(/ID:[^\s]+/, "ID:" + meta['sample'])
                                                    meta['id'] = meta['sample']
                                                    [meta,cram, crai]}
    cram_variant_calling.dump(tag:"cram_variant_calling")

    emit:
    cram_variant_calling    = cram_variant_calling
    versions                = versions
    reports                 = reports
}
