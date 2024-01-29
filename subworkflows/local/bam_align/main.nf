//
// DNA and RNA ALIGNMENT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

// SUBWORKFLOWS
// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT   } from '../bam_convert_samtools/main'
// Map input reads to reference genome in DNA
include { FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP               } from '../fastq_align_bwamem_mem2_dragmap/main'
// Map input reads to reference genome in RNA
include { FASTQ_ALIGN_STAR                              } from '../../nf-core/fastq_align_star/main'
// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                      } from '../bam_merge_index_samtools/main'
// Create samplesheets to restart from mapping
include { CHANNEL_ALIGN_CREATE_CSV                      } from '../channel_align_create_csv/main'
// MODULES
// Run FASTQC
include { FASTQC                                        } from '../../../modules/nf-core/fastqc/main'
// TRIM/SPLIT FASTQ Files
include { FASTP                                         } from '../../../modules/nf-core/fastp/main'


workflow BAM_ALIGN {

    take:
    bwa
    bwamem2
    dragmap
    star_index
    gtf
    input_sample

    main:
    reports   = Channel.empty()
    versions  = Channel.empty()

    // Initialise outputs to emit
    bam_mapped_rna   = Channel.empty()
    bam_mapped_dna   = Channel.empty()
    bam_mapped       = Channel.empty()
    cram_mapped      = Channel.empty()

    // Gather index for mapping given the chosen aligner for DNA
    index_alignment = params.aligner == "bwa-mem" ? bwa :
        params.aligner == "bwa-mem2" ? bwamem2 :
        dragmap
    if (params.step == 'mapping') {

        // Figure out if input is bam or fastq
        input_sample_type = input_sample.branch{
            bam:   it[0].data_type == "bam"
            fastq: it[0].data_type == "fastq"
        }

        // convert any bam input to fastq
        // fasta are not needed when converting bam to fastq -> [ id:"fasta" ], []
        // No need for fasta.fai -> []
        interleave_input = false // Currently don't allow interleaved input
        CONVERT_FASTQ_INPUT(
            input_sample_type.bam,
            [ [ id:"fasta" ], [] ], // fasta
            [ [ id:'null' ], [] ],  // fasta_fai
            interleave_input)

        // Gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        input_fastq = input_sample_type.fastq.mix(CONVERT_FASTQ_INPUT.out.reads)



        // STEP 1.B: QC
        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(input_fastq)

            reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
            versions = versions.mix(FASTQC.out.versions.first())
        }



        //  STEP 1.C: Trimming and/or splitting
        if (params.trim_fastq || params.split_fastq > 0) {

            save_trimmed_fail = false
            save_merged = false
            FASTP(
                input_fastq,
                [], // we are not using any adapter fastas at the moment
                save_trimmed_fail,
                save_merged
            )

            reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
            reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })

            if (params.split_fastq) {
                reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                    read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                    [ meta + [ size:read_files.size() ], read_files ]
                }.transpose()
            } else reads_for_alignment = FASTP.out.reads

            versions = versions.mix(FASTP.out.versions)

        } else {
            reads_for_alignment = input_fastq
        }

        //  STEP 1.D: MAPPING READS TO REFERENCE GENOME
        // Generate mapped reads channel for alignment
        // reads will be sorted
        reads_for_alignment = reads_for_alignment.map{ meta, reads ->
            // Update meta.id to meta.sample no multiple lanes or splitted fastqs
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
            }
        // Separate DNA from RNA samples, DNA samples will be aligned with bwa, and RNA samples with star
        reads_for_alignment_status = reads_for_alignment.branch{
                dna: it[0].status < 2
                rna: it[0].status == 2
            }

        //  STEP 1.D.1: DNA mapping with BWA
        sort_bam = true
        FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP(reads_for_alignment_status.dna, index_alignment, sort_bam)

        // Grouping the bams from the same samples not to stall the workflow
        bam_mapped_dna = FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP.out.bam.map{ meta, bam ->

            // Update meta.id to be meta.sample, ditching sample-lane that is not needed anymore
            // Update meta.data_type
            // Remove no longer necessary fields:
            //   read_group: Now in the BAM header
            //    num_lanes: only needed for mapping
            //         size: only needed for mapping

            // Use groupKey to make sure that the correct group can advance as soon as it is complete
            // and not stall the workflow until all reads from all channels are mapped
            [ groupKey( meta - meta.subMap('num_lanes', 'read_group', 'size') + [ data_type:'bam', id:meta.sample ], (meta.num_lanes ?: 1) * (meta.size ?: 1)), bam ]
        }.groupTuple()
        bam_mapped_dna.dump(tag:"bam_mapped_dna")
        reads_for_alignment_status.rna.dump(tag:"reads_for_alignment_status.rna")

        // RNA will be aligned with STAR
        // Run STAR
        FASTQ_ALIGN_STAR (
            reads_for_alignment_status.rna,
            star_index,
            gtf,
            params.star_ignore_sjdbgtf,
            params.seq_platform ? params.seq_platform : [],
            params.seq_center ? params.seq_center : [],
            [ [ id:"fasta" ], [] ] // fasta
        )
        // Grouping the bams from the same samples not to stall the workflow
        bam_mapped_rna = FASTQ_ALIGN_STAR.out.bam.map{ meta, bam ->

            // Update meta.id to be meta.sample, ditching sample-lane that is not needed anymore
            // Update meta.data_type
            // Remove no longer necessary fields:
            //   read_group: Now in the BAM header
            //    num_lanes: only needed for mapping
            //         size: only needed for mapping

            // Use groupKey to make sure that the correct group can advance as soon as it is complete
            // and not stall the workflow until all reads from all channels are mapped
            [ groupKey( meta - meta.subMap('num_lanes', 'read_group', 'size', 'lane') + [ data_type:'bam', id:meta.sample ], (meta.num_lanes ?: 1) * (meta.size ?: 1)), bam ]
        }.groupTuple()
        bam_mapped_rna.dump(tag:"bam_mapped_rna")
        // Gather QC reports
        reports           = reports.mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]}.ifEmpty([]))
        reports           = reports.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
        versions          = versions.mix(FASTQ_ALIGN_STAR.out.versions)

        // mix dna and rna in one channel
        bam_mapped = bam_mapped_dna.mix(bam_mapped_rna)

        // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
        // Except if and only if skipping markduplicates or saving mapped bams
        if (params.save_mapped || (params.skip_tools && params.skip_tools.split(',').contains('markduplicates'))) {

            // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)
            cram_mapped = BAM_TO_CRAM_MAPPING.out.alignment_index

            // Create CSV to restart from this step
            params.save_output_as_bam ? CHANNEL_ALIGN_CREATE_CSV(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai) : CHANNEL_ALIGN_CREATE_CSV(BAM_TO_CRAM_MAPPING.out.alignment_index)

            // Gather used softwares versions
            versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
            versions = versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
        }

        // Gather used softwares versions
        versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)
        versions = versions.mix(FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP.out.versions)
        versions = versions.mix(FASTQ_ALIGN_STAR.out.versions)

    }

    emit:
    // TODO: do I need to output RNA and DNA separately or cam I directly use bam_mapped but separating them?
    bam_mapped_rna   = bam_mapped_rna  // second pass with RG tags
    bam_mapped_dna   = bam_mapped_dna  // second pass with RG tags
    bam_mapped       = bam_mapped      // for preprocessing
    cram_mapped      = cram_mapped     // for preprocessing
    reports          = reports
    versions         = versions
}
