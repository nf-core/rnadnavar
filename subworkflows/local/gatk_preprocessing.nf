//
// GATK pre-processing best practices
//
// samtools
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM               } from '../../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM_RECAL         } from '../../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM               } from '../../modules/nf-core/modules/samtools/convert/main'
// Mark Duplicates (+QC)
include { MARKDUPLICATES                                       } from '../nf-core/gatk4/markduplicates/main'
include { MARKDUPLICATES_CSV                                   } from '../local/markduplicates_csv'
include { MARKDUPLICATES_SPARK                                 } from '../nf-core/gatk4/markduplicates_spark/main'
// SplitNCigarReads
include { SPLITNCIGAR                                          } from '../nf-core/splitncigar'        // Splits reads that contain Ns in their cigar string
// Convert to CRAM (+QC)
include { BAM_TO_CRAM                                          } from '../nf-core/bam_to_cram'
include { BAM_TO_CRAM as BAM_TO_CRAM_SNCR                      } from '../nf-core/bam_to_cram'
// QC on CRAM
include { CRAM_QC                                              } from '../nf-core/cram_qc'
// Create recalibration tables
include { PREPARE_RECALIBRATION                                } from '../nf-core/gatk4/prepare_recalibration/main'
include { PREPARE_RECALIBRATION_CSV                            } from '../local/prepare_recalibration_csv'
// Create recalibration tables SPARK
include { PREPARE_RECALIBRATION_SPARK                          } from '../nf-core/gatk4/prepare_recalibration_spark/main'
// Create recalibrated cram files to use for variant calling (+QC)
include { RECALIBRATE                                          } from '../nf-core/gatk4/recalibrate/main'
include { RECALIBRATE_SPARK                                    } from '../nf-core/gatk4/recalibrate_spark/main'
include { RECALIBRATE_CSV                                      } from '../local/recalibrate_csv'




workflow GATK_PREPROCESSING {
    take:
        step                          // Mandatory, step to start with
        ch_bam_mapped_dna             // channel: [mandatory] ch_bam_mapped_dna
        ch_bam_mapped_rna             // channel: [mandatory] ch_bam_mapped_rna
        skip_tools                    // channel: [mandatory] skip_tools
        use_gatk_spark                // channel: [mandatory] use_gatk_spark
        save_output_as_bam            // channel: [mandatory] save_output_as_bam
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        dict
        germline_resource             // channel: [optional]  germline_resource
        germline_resource_tbi         // channel: [optional]  germline_resource_tbi
        intervals                     // channel: [mandatory] intervals/target regions
        intervals_for_preprocessing                     // channel: [mandatory] intervals/wes
        ch_interval_list_split
        ch_reports
        ch_versions


    main:
    if (step in ['mapping', 'markduplicates']) {

        // 1. SAMTOOLS_CRAMTOBAM ( to speed up computation)
        // 2. Need fasta for cram compression (maybe just using --fasta, because this reference will be used elsewhere)
        ch_cram_no_markduplicates_restart = Channel.empty()
        ch_cram_markduplicates_no_spark   = Channel.empty()
        ch_cram_markduplicates_spark      = Channel.empty()

// STEP 2: markduplicates (+QC) + convert to CRAM

        // ch_bam_for_markduplicates will countain bam mapped with GATK4_MAPPING when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        // ch_bam_for_markduplicates = step == 'mapping'? ch_bam_mapped : ch_input_sample.map{ meta, input, index -> [meta, input] }

        ch_bam_for_markduplicates = Channel.empty()
        ch_input_cram_indexed     = Channel.empty()

        if (step == 'mapping') {
            // put together DNA and RNA for MARKDUPLICATES if coming from alignment
            ch_bam_for_markduplicates = ch_bam_for_markduplicates.mix(
                                        ch_bam_mapped_dna,
                                        ch_bam_mapped_rna)
            }
        else {
        // input was a BAM and there is no need for alignment
            ch_input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{ch_convert}

            ch_bam_for_markduplicates = ch_convert.bam.map{ meta, bam, bai -> [meta, bam]}

            //In case Markduplicates is run convert CRAM files to BAM, because the tool only runs on BAM files. MD_SPARK does run on CRAM but is a lot slower
            if (!(skip_tools && skip_tools.split(',').contains('markduplicates'))){

                SAMTOOLS_CRAMTOBAM(ch_convert.cram, fasta, fasta_fai)
                ch_versions = ch_versions.mix(SAMTOOLS_CRAMTOBAM.out.versions)

                ch_bam_for_markduplicates = ch_bam_for_markduplicates.mix(SAMTOOLS_CRAMTOBAM.out.alignment_index.map{ meta, bam, bai -> [meta, bam]})
            } else {
                ch_input_cram_indexed     = ch_convert.cram
            }
        }

        if (skip_tools && skip_tools.split(',').contains('markduplicates')) {
            // ch_bam_indexed will countain bam mapped with GATK4_MAPPING when step is mapping
            // which are then merged and indexed
            // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
            ch_bam_indexed = step == 'mapping' ? MERGE_INDEX_BAM.out.bam_bai : ch_convert.bam

            BAM_TO_CRAM(
                ch_bam_indexed,
                ch_input_cram_indexed,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_no_markduplicates_restart = BAM_TO_CRAM.out.cram_converted

            // Gather QC reports
            ch_reports  = ch_reports.mix(BAM_TO_CRAM.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(BAM_TO_CRAM.out.versions)
        }
        else if (use_gatk_spark && use_gatk_spark.contains('markduplicates')) {
            MARKDUPLICATES_SPARK(
                ch_bam_for_markduplicates,
                dict,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)
            ch_cram_markduplicates_spark = MARKDUPLICATES_SPARK.out.cram

            // Gather QC reports
            ch_reports  = ch_reports.mix(MARKDUPLICATES_SPARK.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MARKDUPLICATES_SPARK.out.versions)
        }
        else {
            MARKDUPLICATES(
                ch_bam_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_markduplicates_no_spark = MARKDUPLICATES.out.cram

            // Gather QC reports
            ch_reports  = ch_reports.mix(MARKDUPLICATES.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)
        }

        // ch_md_cram_for_restart contains either:
        // - crams from markduplicates
        // - crams from markduplicates_spark
        // - crams converted from bam mapped when skipping markduplicates
        ch_md_cram_for_restart = Channel.empty().mix(
        ch_cram_markduplicates_no_spark,
        ch_cram_markduplicates_spark,
        ch_cram_no_markduplicates_restart).map {
            meta, cram, crai ->
            //Make sure correct data types are carried through
            [[
                data_type:  "cram",
                id:         meta.id,
                patient:    meta.patient,
                sample:     meta.sample,
                sex:        meta.sex,
                status:     meta.status
                ],
            cram, crai]
            }

//        // CSV should be written for the file actually out, either CRAM or BAM
//        // Create CSV to restart from this step
        if (!(skip_tools && skip_tools.split(',').contains('markduplicates'))) {
            MARKDUPLICATES_CSV(ch_md_cram_for_restart)
            }
    }

        if (step in ['mapping', 'markduplicates', 'splitncigar']) {
            //
            // Need to separate RNA from DNA to do SPLITNCIGARREADS from GATK
            // Logic to separate DNA from RNA samples, DNA samples will be aligned with bwa, and RNA samples with star
            //
            ch_md_cram_for_restart.branch{
                dna: it[0].status < 2
                rna: it[0].status == 2
            }.set{ch_md_cram_for_splitncigar_status}

            // RNA samples only
            ch_md_cram_for_splitncigar = ch_md_cram_for_splitncigar_status.rna
            ch_md_cram_dna = ch_md_cram_for_splitncigar_status.dna

            // TODO: separate this as a step - at the moment concurrent with markduplicates
            // SUBWORKFLOW: SplitNCigarReads from GATK4 over the intervals
            // Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data).
            //

            ch_splitncigar_bam_bai = Channel.empty()
            SPLITNCIGAR (
                ch_md_cram_for_splitncigar,
                fasta,
                fasta_fai,
                dict,
                ch_interval_list_split
            )
            ch_splitncigar_bam_bai  = SPLITNCIGAR.out.bam_bai
//            ch_versions             = ch_versions.mix(SPLITNCIGAR.out.versions.first().ifEmpty(null))
            ch_input_cram_indexed     = Channel.empty() // TODO introduce in a proper way as below
            // SPLINCIGAR produces BAM as output
            BAM_TO_CRAM_SNCR(
                ch_splitncigar_bam_bai,
                ch_input_cram_indexed,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_splitncigar_restart = BAM_TO_CRAM_SNCR.out.cram_converted

            // Gather QC reports
            ch_reports  = ch_reports.mix(BAM_TO_CRAM_SNCR.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(BAM_TO_CRAM_SNCR.out.versions)

            // join again DNA and RNA to continue pre-processing
            ch_cram_for_recalibration = Channel.empty()

            ch_splitncigar_cram_for_restart = ch_cram_for_recalibration.mix(
                                        ch_md_cram_dna,
                                        ch_cram_splitncigar_restart)

            ch_cram_for_recal = ch_splitncigar_cram_for_restart.map{ meta, cram, crai ->
                        //Make sure correct data types are carried through
                        [[
                            data_type:  "cram",
                            id:         meta.id,
                            patient:    meta.patient,
                            sample:     meta.sample,
                            sex:        meta.sex,
                            status:     meta.status
                            ],
                        cram, crai]
                    }



        }

// STEP 3: Create recalibration tables
    if (step in ['mapping', 'markduplicates', 'splitncigar', 'prepare_recalibration']) {
        // TODO change introducing splitncigar
        // Run if starting from step "prepare_recalibration"
        if(step == 'prepare_recalibration'){

            //Support if starting from BAM or CRAM files
            ch_input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{ch_convert}

            //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            SAMTOOLS_BAMTOCRAM(ch_convert.bam, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions)

            ch_cram_for_prepare_recalibration = Channel.empty().mix(SAMTOOLS_BAMTOCRAM.out.alignment_index, ch_convert.cram)

            ch_cram_for_recal = SAMTOOLS_BAMTOCRAM.out.alignment_index

        } else {

            // ch_cram_for_prepare_recalibration contains either:
            // - crams from markduplicates
            // - crams from markduplicates_spark
            // - crams converted from bam mapped when skipping markduplicates
            // - input cram files, when start from step markduplicates
            //ch_cram_for_recal.view() //contains md.cram.crai or sncr.cram.crai
            ch_cram_for_prepare_recalibration = Channel.empty().mix(ch_cram_for_recal, ch_input_cram_indexed)
        }

        if (!(skip_tools && skip_tools.split(',').contains('baserecalibrator'))) {
            ch_table_bqsr_no_spark = Channel.empty()
            ch_table_bqsr_spark    = Channel.empty()

            if (use_gatk_spark && use_gatk_spark.contains('baserecalibrator')) {
                PREPARE_RECALIBRATION_SPARK(
                    ch_cram_for_prepare_recalibration,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals,
                    germline_resource,
                    germline_resource_tbi)

                    ch_table_bqsr_spark = PREPARE_RECALIBRATION_SPARK.out.table_bqsr

                    // Gather used softwares versions
                    ch_versions = ch_versions.mix(PREPARE_RECALIBRATION_SPARK.out.versions)
            } else {

                PREPARE_RECALIBRATION(
                    ch_cram_for_prepare_recalibration,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals,
                    germline_resource,
                    germline_resource_tbi)

                    ch_table_bqsr_no_spark = PREPARE_RECALIBRATION.out.table_bqsr

                    // Gather used softwares versions
                    ch_versions = ch_versions.mix(PREPARE_RECALIBRATION.out.versions)
            }

            // ch_table_bqsr contains either:
            // - bqsr table from baserecalibrator
            // - bqsr table from baserecalibrator_spark
            ch_table_bqsr = Channel.empty().mix(
                ch_table_bqsr_no_spark,
                ch_table_bqsr_spark)

            ch_reports  = ch_reports.mix(ch_table_bqsr.collect{ meta, table -> table})

            ch_cram_applybqsr = ch_cram_for_prepare_recalibration.join(ch_table_bqsr)

            // Create CSV to restart from this step
            PREPARE_RECALIBRATION_CSV(ch_cram_for_recal.join(ch_table_bqsr), skip_tools)
        }
    }
// STEP 4: RECALIBRATING
    if (step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        // Run if starting from step "prepare_recalibration"
        if(step == 'recalibrate'){
            //Support if starting from BAM or CRAM files
            ch_input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }.set{ch_convert}

            //If BAM file, split up table and mapped file to convert BAM to CRAM
            ch_bam_table = ch_convert.bam.map{ meta, bam, bai, table -> [meta, table]}
            ch_bam_bam   = ch_convert.bam.map{ meta, bam, bai, table -> [meta, bam, bai]}

            //BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            SAMTOOLS_BAMTOCRAM(ch_bam_bam, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions)

            ch_cram_applybqsr = Channel.empty().mix(
                                    SAMTOOLS_BAMTOCRAM.out.alignment_index.join(ch_bam_table),
                                    ch_convert.cram) // Join together converted cram with input tables
        }

        if (!(skip_tools && skip_tools.split(',').contains('baserecalibrator'))) {
            ch_cram_variant_calling_no_spark = Channel.empty()
            ch_cram_variant_calling_spark    = Channel.empty()

            if (use_gatk_spark && use_gatk_spark.contains('baserecalibrator')) {
                RECALIBRATE_SPARK(
                    ch_cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals)

                ch_cram_variant_calling_spark = RECALIBRATE_SPARK.out.cram

                // Gather used softwares versions
                ch_versions = ch_versions.mix(RECALIBRATE_SPARK.out.versions)

            } else {
                RECALIBRATE(
                    ch_cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals)

                ch_cram_variant_calling_no_spark = RECALIBRATE.out.cram

                // Gather used softwares versions
                ch_versions = ch_versions.mix(RECALIBRATE.out.versions)
            }
            ch_cram_variant_calling = Channel.empty().mix(
                ch_cram_variant_calling_no_spark,
                ch_cram_variant_calling_spark)

            CRAM_QC(
                ch_cram_variant_calling,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            // Gather QC reports
            ch_reports  = ch_reports.mix(CRAM_QC.out.qc.collect{meta, report -> report})

            // Gather used softwares versions
            ch_versions = ch_versions.mix(CRAM_QC.out.versions)

            //If save_output_as_bam, then convert CRAM files to BAM
            SAMTOOLS_CRAMTOBAM_RECAL(ch_cram_variant_calling, fasta, fasta_fai)
            ch_versions = ch_versions.mix(SAMTOOLS_CRAMTOBAM_RECAL.out.versions)

            // CSV should be written for the file actually out out, either CRAM or BAM
            csv_recalibration = Channel.empty()
            csv_recalibration = save_output_as_bam ?  SAMTOOLS_CRAMTOBAM_RECAL.out.alignment_index : ch_cram_variant_calling

            // Create CSV to restart from this step
            RECALIBRATE_CSV(csv_recalibration)


        } else if (step == 'recalibrate'){
            // ch_cram_variant_calling contains either:
            // - input bams converted to crams, if started from step recal + skip BQSR
            // - input crams if started from step recal + skip BQSR
            ch_cram_variant_calling = Channel.empty().mix(SAMTOOLS_BAMTOCRAM.out.alignment_index,
                                                        ch_convert.cram.map{ meta, cram, crai, table -> [meta, cram, crai]})
        } else {
            // ch_cram_variant_calling contains either:
            // - crams from markduplicates = ch_cram_for_prepare_recalibration if skip BQSR but not started from step recalibration
            ch_cram_variant_calling = Channel.empty().mix(ch_cram_for_prepare_recalibration)
        }
    }
// SEPARATE ABOVE AS A SUBWORKFLOW - GATK_PREPROCESSING
    
    
    emit:
    ch_cram_variant_calling

    versions    = ch_versions
    ch_reports = ch_reports
}
    
    