//
// GATK pre-processing best practices
//
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM               } from '../../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM_RECAL         } from '../../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_CRAMTOBAM_SNCR          } from '../../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM               } from '../../modules/nf-core/modules/samtools/convert/main'
include { MARKDUPLICATES                                       } from '../nf-core/gatk4/markduplicates/main'
include { MARKDUPLICATES_CSV                                   } from '../local/markduplicates_csv'
include { SPLITNCIGAR                                          } from '../nf-core/splitncigar'        // Splits reads that contain Ns in their cigar string
include { BAM_TO_CRAM                                          } from '../nf-core/bam_to_cram'
include { BAM_TO_CRAM as BAM_TO_CRAM_SNCR                      } from '../nf-core/bam_to_cram'
include { CRAM_QC                                              } from '../nf-core/cram_qc'
include { PREPARE_RECALIBRATION                                } from '../nf-core/gatk4/prepare_recalibration/main'
include { PREPARE_RECALIBRATION_CSV                            } from '../local/prepare_recalibration_csv'
include { RECALIBRATE                                          } from '../nf-core/gatk4/recalibrate/main'
include { RECALIBRATE_CSV                                      } from '../local/recalibrate_csv'

workflow GATK_PREPROCESSING {
    take:
        step                          // Mandatory, step to start with
        ch_bam_mapped                 // channel: [mandatory] ch_bam_mapped
        skip_tools                    // channel: [mandatory] skip_tools
        save_output_as_bam            // channel: [mandatory] save_output_as_bam
        fasta                         // channel: [mandatory] fasta
        fasta_fai                     // channel: [mandatory] fasta_fai
        dict
        germline_resource             // channel: [optional]  germline_resource
        germline_resource_tbi         // channel: [optional]  germline_resource_tbi
        intervals                     // channel: [mandatory] intervals/target regions
        intervals_for_preprocessing   // channel: [mandatory] intervals/wes
        ch_interval_list_split
        ch_input_sample


    main:
        ch_reports   = Channel.empty()
        ch_versions  = Channel.empty()

        // Select inputs for makduplicates/recalibration
        ch_bam_for_markduplicates         = Channel.empty()
        ch_input_cram_indexed             = Channel.empty()
        ch_cram_no_markduplicates_restart = Channel.empty()
        ch_cram_markduplicates            = Channel.empty()
        ch_cram_variant_calling           = Channel.empty()
        // input from mapping
        if (step == 'mapping' | !ch_input_sample) {
            ch_bam_for_markduplicates = ch_bam_mapped
            ch_input_sample = ch_bam_mapped
        } else {
                // input from samplesheet was a BAM and there is no need for alignment
                ch_bam_mapped.branch{
                    bam:  it[0].data_type == "bam"
                    cram: it[0].data_type == "cram"
                }.set{ch_convert}
                ch_bam_for_markduplicates = ch_convert.bam.map{ meta, bam, bai -> [meta, bam]}
                // If CRAM files, convert to BAM, because the tool only runs on BAM files.
                if (!(skip_tools && skip_tools.split(',').contains('markduplicates'))){
                    // SAMTOOLS_CRAMTOBAM ( to speed up computation)
                    SAMTOOLS_CRAMTOBAM(ch_convert.cram, fasta, fasta_fai)
                    ch_versions = ch_versions.mix(SAMTOOLS_CRAMTOBAM.out.versions)
                    ch_bam_for_markduplicates = ch_bam_for_markduplicates.mix(SAMTOOLS_CRAMTOBAM.out.alignment_index
                                                                         .map{ meta, bam, bai -> [meta, bam]})
                } else {
                    ch_cram_no_markduplicates_restart  = ch_convert.cram
                    ch_cram_markduplicates             = Channel.empty()
                    // ch_bam_for_markduplicates will countain bam mapped with GATK4_MAPPING when step is mapping
                    // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration

                    ch_bam_for_markduplicates = ch_convert.bam
                    ch_input_cram_indexed     = Channel.empty()
                }
        }

        // STEP 1: mark duplicates
        if (step in ['mapping', 'markduplicates'] ) {
            // NO markduplicates, just convert BAM to CRAM
            if (skip_tools && skip_tools.split(',').contains('markduplicates')) {
                // ch_bam_indexed will countain bam mapped with GATK4_MAPPING when step is mapping
                // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
                ch_bam_indexed = step == 'mapping' ? ch_bam_mapped : ch_convert.bam
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
            // MARKDUPLICATES
            else {
//                ch_bam_for_markduplicates = ch_bam_for_markduplicates.map{meta, bams ->[ meta, bams.sort()]}
                MARKDUPLICATES(
                    ch_bam_for_markduplicates,
                    fasta,
                    fasta_fai,
                    intervals_for_preprocessing)
                ch_cram_markduplicates = MARKDUPLICATES.out.cram
                // Gather QC reports
                ch_reports  = ch_reports.mix(MARKDUPLICATES.out.qc.collect{meta, report -> report})
                // Gather used softwares versions
                ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions)
                }
            // ch_md_cram_for_restart contains either:
            // - crams from markduplicates
            // - crams converted from bam mapped when skipping markduplicates
            ch_md_cram_for_restart = Channel.empty()
                                            .mix(
                                                ch_cram_markduplicates,
                                                ch_cram_no_markduplicates_restart)
                                            .map {
                                                meta, cram, crai ->
                                                //Make sure correct data types are carried through
                                                [[
                                                    data_type:  "cram",
                                                    id:         meta.id,
                                                    patient:    meta.patient,
                                                    sample:     meta.sample,
                                                    status:     meta.status,
                                                    lane:       meta.lane
                                                    ],
                                                cram, crai]
                                                }
            // CSV should be written for the file actually out, either CRAM or BAM
            // Create CSV to restart from this step
            if (!(skip_tools && skip_tools.split(',').contains('markduplicates'))) {
                MARKDUPLICATES_CSV(ch_md_cram_for_restart)
            }
        }
        // STEP 1b: SplitNCigarReads for RNA
        if (step in ['mapping', 'markduplicates', 'splitncigar']) {
            if (step == 'splitncigar') {
            ch_md_cram_for_restart = ch_bam_mapped
            }

            // Separate RNA from DNA to do SPLITNCIGARREADS from GATK
            ch_md_cram_for_restart.branch{
                dna: it[0].status < 2
                rna: it[0].status == 2
            }.set{ch_md_for_splitncigar}
            // RNA samples only
            ch_for_splitncigar = ch_md_for_splitncigar.rna
            ch_md_cram_dna = ch_md_for_splitncigar.dna
            // If CRAM files, convert to BAM, because the tool only runs on BAM files.
            ch_for_splitncigar.branch{
                    bam:  it[0].data_type == "bam"
                    cram: it[0].data_type == "cram"
                }.set{ch_for_splitncigar_input}
            SAMTOOLS_CRAMTOBAM_SNCR(ch_for_splitncigar_input.cram, fasta, fasta_fai)
            ch_md_cram_for_splitncigar = ch_for_splitncigar_input.bam.mix(SAMTOOLS_CRAMTOBAM_SNCR.out.alignment_index)
            SPLITNCIGAR (
                ch_md_cram_for_splitncigar,
                fasta,
                fasta_fai,
                dict,
                intervals_for_preprocessing
            )
            ch_splitncigar_bam_bai  = SPLITNCIGAR.out.bam_bai
            ch_versions             = ch_versions.mix(SPLITNCIGAR.out.versions)
            // empty channel as BAM_TO_CRAM needs it as input
            ch_input_cram_indexed     = Channel.empty()
            // SPLINCIGAR BAM to CRAM
            BAM_TO_CRAM_SNCR(
                ch_splitncigar_bam_bai,
                ch_input_cram_indexed,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            ch_cram_splitncigar = BAM_TO_CRAM_SNCR.out.cram_converted
            // Gather QC reports
            ch_reports  = ch_reports.mix(BAM_TO_CRAM_SNCR.out.qc.collect{meta, report -> report})
            // Gather used softwares versions
            ch_versions = ch_versions.mix(BAM_TO_CRAM_SNCR.out.versions)
            // join again DNA and RNA to continue pre-processing
            ch_cram_for_recalibration = Channel.empty()
            ch_splitncigar_cram_for_restart = ch_cram_for_recalibration.mix(
                                        ch_md_cram_dna,
                                        ch_cram_splitncigar)
            ch_cram_for_recal = ch_splitncigar_cram_for_restart.map{ meta, cram, crai ->
                        [[
                            data_type:  "cram",
                            id:         meta.id,
                            patient:    meta.patient,
                            sample:     meta.sample,
                            status:     meta.status
                            ],
                        cram, crai]
                    }
        }
        // STEP 2: Create recalibration tables
        if (step in ['mapping', 'markduplicates', 'splitncigar', 'prepare_recalibration'] ) {
            // Run if starting from step "prepare_recalibration" or from "splitncigar" but skipping splitncigar (not much sense but just in case)
            if (step  == 'prepare_recalibration' || (step == "splitncigar" && skip_tools.split(',').contains('splitncigar'))){
                // Known issue: if you try to start the pipeline with a csv with a file in the table column
                // but want to start from prepare_recalibration (re-do the table) then it will throw an error
                // because the input object already has a table - this is actually bad practice for the pipeline so should not be used.
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

            }
            else {
                // ch_cram_for_prepare_recalibration contains either:
                // - crams converted from bam mapped when skipping markduplicates
                // - input cram files, when start from step markduplicates
                ch_cram_for_prepare_recalibration = Channel.empty().mix(ch_cram_for_recal, ch_input_cram_indexed)
            }
            // BASERECALIBRATOR
            if (!(skip_tools &&  skip_tools.split(',').contains('baserecalibrator'))) {
                ch_table_bqsr = Channel.empty()

                ch_cram_for_prepare_recalibration.dump(tag:"[STEP2_GATKPREPROCESSING] cram_input_for_recal")
                PREPARE_RECALIBRATION(
                    ch_cram_for_prepare_recalibration,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals,
                    germline_resource,
                    germline_resource_tbi)

                ch_table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr
                // Gather used softwares versions
                ch_versions = ch_versions.mix(PREPARE_RECALIBRATION.out.versions)


                ch_reports  = ch_reports.mix(ch_table_bqsr.collect{ meta, table -> table})
                ch_cram_applybqsr = ch_cram_for_prepare_recalibration.join(ch_table_bqsr)
                // Create CSV to restart from this step
                PREPARE_RECALIBRATION_CSV(ch_cram_for_recal.join(ch_table_bqsr), skip_tools)
            }
        }
        // STEP 3: RECALIBRATING
        if (step in ['mapping', 'markduplicates', 'splitncigar', 'prepare_recalibration', 'recalibrate'] ) {
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

                ch_bam_table.dump(tag:"ch_bam_table")
                ch_bam_bam.dump(tag:"ch_bam_bam")
//                BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
                SAMTOOLS_BAMTOCRAM(ch_bam_bam, fasta, fasta_fai)
                ch_versions = ch_versions.mix(SAMTOOLS_BAMTOCRAM.out.versions)

                ch_cram_applybqsr = Channel.empty().mix(
                                        SAMTOOLS_BAMTOCRAM.out.alignment_index.join(ch_bam_table),
                                        ch_convert.cram) // Join together converted cram with input tables
            }
        if (!(skip_tools && skip_tools.split(',').contains('baserecalibrator'))) {
            // RECALIBRATION
//            ch_cram_applybqsr.dump(tag:"[STEP2_GATKPREPROCESSING] PREPARE_RECALIBRATION")

            RECALIBRATE(
                ch_cram_applybqsr,
                dict,
                fasta,
                fasta_fai,
                intervals)

            ch_cram_variant_calling = RECALIBRATE.out.cram
            ch_versions = ch_versions.mix(RECALIBRATE.out.versions)

            // QC for resulting CRAM(s)
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
            RECALIBRATE_CSV(csv_recalibration.transpose())

            } else if (step == 'recalibrate'){
                // ch_cram_variant_calling contains either:
                // - input bams converted to crams, if started from step recal + skip BQSR
                // - input crams if started from step recal + skip BQSR
                ch_cram_variant_calling = Channel.empty().mix(SAMTOOLS_BAMTOCRAM.out.alignment_index,
                                                            ch_convert.cram.map{ meta, cram, crai, table -> [meta, cram, crai]})
            } else {
                // ch_cram_variant_calling contains either:
                // - crams from markduplicates = ch_cram_for_prepare_recalibration if skip BQSR but not started from step recalibration
                ch_cram_variant_calling = ch_cram_for_prepare_recalibration
            }
        }


    emit:
        ch_cram_variant_calling = ch_cram_variant_calling
        versions    = ch_versions
        ch_reports = ch_reports
}
    
    