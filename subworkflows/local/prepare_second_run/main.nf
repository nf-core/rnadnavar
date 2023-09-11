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
// Realignment with HISAT2
include { FASTQ_ALIGN_HISAT2                            } from '../../nf-core/fastq_align_hisat2/main'


workflow BAM_EXTRACT_READS_HISAT2_ALIGN {
    take:
        input_sample
        maf_with_candidates              // MAf with candidate regions to extract
        reads_to_realign                 // CRAM/BAM to extract reads from
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
        'prepare_recalibration', 'recalibrate', 'variant_calling', 'normalize', 'consensus', 
        'second_run'] && !(params.skip_tools && params.skip_tools.split(",").contains("second_run"))) {
	        if (params.step == 'second_run') {
	            input_elements_status = input_sample.branch{
	                                    norealign: it[0].status == 1
				                        realign:   it[0].status == 2 || it[0].status == 0
	            }
				input_elements_status.norealign.dump(tag:"input_elements_status.norealign")
	            dna_mafs = input_elements_status.norealign.map{meta, vcf -> [meta + [ consensus: meta.variantcaller ==~ /(\S+)?(?i)consensus(\S+)/ ], vcf]}
	            dna_mafs_consensus = dna_mafs.branch{
	                                    isconsensus: it[0].consensus == true
				                        noconsensus: it[0].consensus == false
	                                     }
	            dna_consensus_maf = dna_mafs_consensus.isconsensus
                dna_varcall_mafs  = dna_mafs_consensus.noconsensus
	            // [meta, cram, crai, maf] (RNA dna NORMAL)
	            cram_to_realign = input_elements_status.realign
	            // TODO convert to CRAM or/and MERGE if necessary
	            // TODO Merge alignments if applicable
//	            MERGE_ALIGN(previous_alignment.map{meta, cram, crai, maf -> [meta, cram]})
	        }  else {
	            // TODO (remember to change id to _realign)
	            cram_to_realign = Channel.empty()
	        }
	        // Get candidate regions
	        // Add files to meta to keep them for next processes
	        maf_to_bed = cram_to_realign.map{meta, cram, crai, maf -> [meta + [cram_file:cram, crai_file:crai, maf_file:maf], maf]}
	        maf_to_bed.dump(tag:"maf_to_bed")
	        MAF2BED(maf_to_bed)
	        // Extract read names with regions from bed
            cram_to_extract = MAF2BED.out.bed.map{meta, bed -> [meta, meta.cram_file, meta.crai_file, bed]}
            cram_to_extract.dump(tag:"cram_to_extract")
	        SAMTOOLS_EXTRACT_READ_IDS(cram_to_extract)
	        // Extract reads
	        cram_to_filter = SAMTOOLS_EXTRACT_READ_IDS.out.read_ids.map{meta, readsid -> [meta, meta.cram_file, readsid]}
	        cram_to_filter.dump(tag:"cram_to_filter")
            PICARD_FILTERSAMREADS(cram_to_filter, 'includeReadList') // bam -> filtered_bam
            // Conver to FQ
            bam_to_fq = PICARD_FILTERSAMREADS.out.bam.join(PICARD_FILTERSAMREADS.out.bai)
            bam_to_fq.dump(tag:"bam_to_fq")
            fasta.dump(tag:"fasta")
            fasta_fai.dump(tag:"fasta_fai")
            interleave_input = false // Currently don't allow interleaved input
            CONVERT_FASTQ_INPUT(
                                bam_to_fq,
                                fasta.map{it -> [ [ id:"fasta" ], it ]}, // fasta
                                fasta_fai.map{it -> [ [ id:"fasta_fai" ], it ]},  // fasta_fai
                                interleave_input
                                )
			// Align with HISAT2
            reads_for_realignment = CONVERT_FASTQ_INPUT.out.reads
            reads_for_realignment.dump(tag:"reads_for_realignment")
            // TODO: add single_end to input check
            FASTQ_ALIGN_HISAT2(
                                reads_for_realignment.map{meta, reads -> [meta + [single_end:false], reads]},
                                hisat2_index,
                                splicesites,
                                fasta.map{it -> [ [ id:"fasta" ], it ]}
                                )
            // Mix with index
            bam_mapped = FASTQ_ALIGN_HISAT2.out.bam.mix(FASTQ_ALIGN_HISAT2.out.bai)
    }

    emit:
    bam_mapped         = bam_mapped
    dna_consensus_maf  = dna_consensus_maf
    dna_varcall_mafs   = dna_varcall_mafs
    versions           = versions                                                         // channel: [ versions.yml ]








}