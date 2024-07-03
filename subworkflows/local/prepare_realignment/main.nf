//
// PREPARE SECOND RUN: extract reads from candidate regions for re-alignment (RNA and DNA normal only)
//
include { MAF2BED                                      } from '../../../modules/local/maf2bed/main'
// Extract read ids for selected regions
include { SAMTOOLS_EXTRACT_READ_IDS                    } from '../../../modules/local/extract_reads_id/main'
// Filter bam for selected regions
include { PICARD_FILTERSAMREADS                        } from '../../../modules/nf-core/picard/filtersamreads/main'
// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT  } from '../bam_convert_samtools/main'
// Conver CRAM to BAM (picard/filtersamreads can't work with cram)
include { SAMTOOLS_CONVERT as CONVERT_CRAM2BAM         } from '../../../modules/nf-core/samtools/convert/main'
// Realignment with HISAT2
include { FASTQ_ALIGN_HISAT2                            } from '../../nf-core/fastq_align_hisat2/main'


workflow BAM_EXTRACT_READS_HISAT2_ALIGN {
    take:
        input_sample
        maf_with_candidates              // MAf with candidate regions to extract [meta, maf]
        reads_to_realign                 // CRAM/BAM to extract reads from [meta, cram, crai]
        fasta
        fasta_fai
        dict
        hisat2_index
        splicesites
        dna_consensus_maf
        dna_varcall_mafs

    main:
        versions   = Channel.empty()
        bam_mapped = Channel.empty()

        if (params.step in ['mapping', 'markduplicates', 'splitncigar',
        'prepare_recalibration', 'recalibrate', 'variant_calling', 'norm', 'consensus',
        'realignment'] && !(params.skip_tools && params.skip_tools.split(",").contains("realignment"))) {
            if (params.step == 'realignment') {
                input_elements_status = input_sample.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
                }  // [meta, maf]
                input_elements_status.norealign.dump(tag:"input_elements_status.norealign")
                dna_mafs = input_elements_status.norealign.map{meta, vcf -> [meta + [ consensus: meta.variantcaller ==~ /(\S+)?(?i)consensus(\S+)/ ], vcf, meta.variantcaller]}
                dna_mafs_consensus = dna_mafs.branch{
                                        isconsensus: it[0].consensus == true
                                        noconsensus: it[0].consensus == false
                                        }
                dna_consensus_maf = dna_mafs_consensus.isconsensus
                dna_varcall_mafs  = dna_mafs_consensus.noconsensus.map{meta, maf, vc -> [meta, [maf], [vc]]}
                dna_mafs.dump(tag:"dna_mafs")
                dna_consensus_maf.dump(tag:"dna_consensus_maf")
                maf_with_candidates.dump(tag:"maf_with_candidates")
                // [meta, cram, crai, maf] (RNA dna NORMAL)
                cram_to_realign = input_elements_status.realign  // [meta, cram, crai, maf]
                cram_to_realign.dump(tag:"cram_to_realign")
                // TODO convert to CRAM or/and MERGE if necessary
                // TODO Merge alignments if applicable
//	            MERGE_ALIGN(previous_alignment.map{meta, cram, crai, maf -> [meta, cram]})
            }  else {
                reads_to_realign.dump(tag:'reads_to_realign0')
                // Filter DNA tumour samples, only DNA matched normal and RNA will be processed for the realignment
                reads_to_realign_branch = reads_to_realign.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
                                        }  // [meta, reads]
                maf_with_candidates_branch = maf_with_candidates.branch{
                                        norealign: it[0].status == 1
                                        realign:   it[0].status == 2 || it[0].status == 0
                                        }
                maf_with_candidates_branch.realign.dump(tag:'maf_with_candidates0')
                // map and change ids with _realign tag to differentiate from first alignment
                maf_with_candidates_tumor = maf_with_candidates_branch.realign.map{
                                        meta, maf ->
                                        [[
                                        patient: meta.patient,
                                        sample:  meta.id.split('_vs_')[0].split('_with_')[0] + "_realign",
                                        status:  meta.status,
                                        id:      meta.id.split('_vs_')[0].split('_with_')[0] + "_realign"
                                        ], maf]
                                        }
                maf_with_candidates_normal = maf_with_candidates_branch.realign.map{
                                        meta, maf ->
                                        [[
                                        patient: meta.patient,
                                        sample:  meta.id.split('_vs_')[1].split('_with_')[0] + "_realign",
                                        status:  0,
                                        id:      meta.id.split('_vs_')[1].split('_with_')[0] + "_realign"
                                        ], maf]
                                        }
                maf_with_candidates_to_realign = maf_with_candidates_tumor.mix(maf_with_candidates_normal)
                reads_to_realign_and_join = reads_to_realign_branch.realign.map{meta, cram, crai ->
                                        [[patient: meta.patient,
                                        sample:  meta.id + "_realign",
                                        status:  meta.status,
                                        id:      meta.id + "_realign"], cram, crai]
                                        }
                reads_to_realign_and_join.dump(tag:'reads_to_realign1')
                maf_with_candidates_to_realign.dump(tag:'maf_with_candidates_to_realign')
                cram_to_realign = reads_to_realign_and_join.join(maf_with_candidates_to_realign)
                cram_to_realign.dump(tag:"cram_to_realign")
            }
            // Get candidate regions
            // Add files to meta to keep them for next processes
            maf_to_bed = cram_to_realign.map{meta, cram, crai, maf -> [meta + [cram_file:cram, crai_file:crai, maf_file:maf], maf]}
            MAF2BED(maf_to_bed)
            // Extract read names with regions from bed
            cram_to_extract = MAF2BED.out.bed.map{meta, bed -> [meta, meta.cram_file, meta.crai_file, bed]}
            SAMTOOLS_EXTRACT_READ_IDS(cram_to_extract)
            // Extract reads
            cram_to_convert = SAMTOOLS_EXTRACT_READ_IDS.out.read_ids.map{meta, readsid -> [meta + [readsid_file:readsid], meta.cram_file, meta.crai_file]}
            // 1) Convert cram 2 bam
            CONVERT_CRAM2BAM(cram_to_convert, fasta, fasta_fai.map{fai -> [[id:"fai"], fai]})
            bam_to_filter = CONVERT_CRAM2BAM.out.bam.map{meta, bam -> [meta, bam, meta.readsid_file]}
            // 2) Apply picard filtersamreads
            PICARD_FILTERSAMREADS(bam_to_filter, fasta,'includeReadList') // bam -> filtered_bam
            // Conver to FQ
            bam_to_fq = PICARD_FILTERSAMREADS.out.bam.join(PICARD_FILTERSAMREADS.out.bai)
            interleave_input = false // Currently don't allow interleaved input
            CONVERT_FASTQ_INPUT(
                                bam_to_fq,
                                fasta,
                                fasta_fai.map{it -> [ [ id:"fasta_fai" ], it ]},
                                interleave_input
                                )
            // Align with HISAT2
            reads_for_realignment = CONVERT_FASTQ_INPUT.out.reads
            hisat2_index.dump(tag:"HISAT2index")
            splicesites.dump(tag:"HISAT2splicesites")
            // Note: single_end in meta always false for this subworkflow TODO: add to samplesheet in future?
            FASTQ_ALIGN_HISAT2(
                                reads_for_realignment.map{meta, reads -> [meta + [single_end:false], reads]},
                                hisat2_index,
                                splicesites,
                                fasta
                                )
            // Mix with index add data type and change id to sample
            bam_mapped = FASTQ_ALIGN_HISAT2.out.bam.join(FASTQ_ALIGN_HISAT2.out.bai).map{meta,bam,bai -> [meta + [ id:meta.sample, data_type:"bam"], bam, bai]}
    }
    bam_mapped = bam_mapped.map{meta, bam, bai -> [meta - meta.subMap('single_end'), bam]}

    emit:
    bam_mapped         = bam_mapped
    dna_consensus_maf  = dna_consensus_maf
    dna_varcall_mafs   = dna_varcall_mafs
    versions           = versions                                                         // channel: [ versions.yml ]

}
