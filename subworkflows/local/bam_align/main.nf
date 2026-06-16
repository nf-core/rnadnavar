//
// DNA and RNA ALIGNMENT
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// SUBWORKFLOWS
// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT   } from '../bam_convert_samtools/main'
// Map input reads to reference genome in DNA
include { FASTQ_ALIGN                                   } from '../fastq_align/main'
// Map input reads to reference genome in RNA
include { FASTQ_ALIGN_STAR                              } from '../../nf-core/fastq_align_star/main'
// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                      } from '../bam_merge_index_samtools/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING       } from '../../../modules/nf-core/samtools/convert/main'
// Create samplesheets to restart from mapping
include { CHANNEL_ALIGN_CREATE_CSV                      } from '../channel_align_create_csv/main'
// MODULES: QC & trimming
include { FASTQC                                        } from '../../../modules/nf-core/fastqc/main'
include { FASTP                                         } from '../../../modules/nf-core/fastp/main'


workflow BAM_ALIGN {

    take:
    bwa
    bwamem2
    dragmap
    star_index
    gtf
    fasta
    fasta_fai
    input_sample

    main:
    reports   = Channel.empty()
    versions  = Channel.empty()
    fasta_with_fai = fasta.combine(fasta_fai).map { meta, fa, fai -> [meta, fa, fai] }

    // Initialize outputs to emit
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
        input_sample_type = input_sample.branch { item ->
            bam:   item[0].data_type == "bam"
            fastq: item[0].data_type == "fastq"
        }
        // QC & TRIM
        interleave_input = false // Currently don't allow interleaved input
        CONVERT_FASTQ_INPUT(
            input_sample_type.bam,
            fasta, // fasta
            fasta_fai,  // fasta_fai
            interleave_input)

        // Gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        input_fastq = input_sample_type.fastq.mix(CONVERT_FASTQ_INPUT.out.reads)
        // QC
        if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(input_fastq)

            reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
            versions = versions.mix(FASTQC.out.versions_fastqc)
        }
        //  Trimming and/or splitting
        if (params.trim_fastq || params.split_fastq > 0) {

            save_trimmed_fail = false
            save_merged = false
            FASTP(
                // we are not using any adapter fastas at the moment, that's why last element is empty
                input_fastq.map { meta, reads -> [meta, reads, []] }, 
                false, // we don't use discard_trimmed_pass at the moment
                save_trimmed_fail,
                save_merged
            )

            reports = reports.mix(FASTP.out.json.collect{ _meta, json -> json })
            reports = reports.mix(FASTP.out.html.collect{ _meta, html -> html })

            if (params.split_fastq) {
                reads_for_alignment = FASTP.out.reads.map { item ->
                    def meta = item[0]
                    def reads = item[1]
                    def read_files = reads.sort(false) { a, b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                    [ meta + [ n_fastq: read_files.size() ], read_files ]
                }.transpose()
            } else reads_for_alignment = FASTP.out.reads

            versions = versions.mix(FASTP.out.versions_fastp)

        } else {
            reads_for_alignment = input_fastq
        }


        //  STEP 1: MAPPING READS TO REFERENCE GENOME
        // Generate mapped reads channel for alignment
        // reads will be sorted
        // First, we must calculate number of lanes for each sample (meta.n_fastq)
        // This is needed to group reads from the same sample together using groupKey to avoid stalling the workflow
        // when reads from different samples are mixed together
        reads_for_alignment.map { item ->
                def meta = item[0]
                def reads = item[1]
                [ meta.subMap('patient', 'sample', 'status'), reads ]
            }
            .groupTuple()
            .map { item ->
                def meta = item[0]
                def reads = item[1]
                meta + [ n_fastq: reads.size() ] // We can drop the FASTQ files now that we know how many there are
            }
            .set { reads_grouping_key }

        reads_for_alignment = reads_for_alignment.map { item ->
            def meta = item[0]
            def reads = item[1]
            // Update meta.id to meta.sample no multiple lanes or splitted fastqs
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
        }
        // Separate DNA from RNA samples, DNA samples will be aligned with bwa, and RNA samples with star
        reads_for_alignment_status = reads_for_alignment.branch { item ->
                dna: item[0].status < 2
                rna: item[0].status == 2
            }

        //  DNA mapping
        sort_bam = true
        FASTQ_ALIGN(reads_for_alignment_status.dna, index_alignment, sort_bam)
        FASTQ_ALIGN.out.bam.dump(tag:'FASTQ_ALIGN.out.bam')
        // Grouping the bams from the same samples not to stall the workflow
        bam_mapped_dna = FASTQ_ALIGN.out.bam
            .combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
            .filter { meta1, _bam, meta2 -> meta1.sample == meta2.sample }
            // Add n_fastq and other variables to meta
            .map { meta1, bam, meta2 ->
                [ meta1 + meta2, bam ]
            }
            // Manipulate meta map to remove old fields and add new ones
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size') + [ data_type: 'bam', id: meta.sample ], bam ]
            }
            // Create groupKey from meta map
            .map { meta, bam ->
                [ groupKey( meta, meta.n_fastq), bam ]
            }
            // Group
            .groupTuple()

        bam_mapped_dna.dump(tag:"bam_mapped_dna")
        reads_for_alignment_status.rna.dump(tag:"reads_for_alignment_status.rna")

        // RNA STAR alignment
        FASTQ_ALIGN_STAR (
            reads_for_alignment_status.rna,
            star_index,
            gtf,
            params.star_ignore_sjdbgtf,
            fasta_with_fai,
            [ [ id:"transcript_fasta" ], [], [] ]
        )
        // Grouping the bams from the same samples not to stall the workflow
        bam_mapped_rna = FASTQ_ALIGN_STAR.out.bam.combine(reads_grouping_key) // Creates a tuple of [ meta, bam, reads_grouping_key ]
            .filter { meta1, _bam, meta2 -> meta1.sample == meta2.sample }
            // Add n_fastq and other variables to meta
            .map { meta1, bam, meta2 ->
                [ meta1 + meta2, bam ]
            }
            // Manipulate meta map to remove old fields and add new ones
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'data_type', 'num_lanes', 'read_group', 'size') + [ data_type: 'bam', id: meta.sample ], bam ]
            }
            // Create groupKey from meta map
            .map { meta, bam ->
                [ groupKey( meta, meta.n_fastq), bam ]
            }
            // Group
            .groupTuple()
        bam_mapped_rna.dump(tag:"bam_mapped_rna")
        // Gather QC reports
        reports           = reports.mix(FASTQ_ALIGN_STAR.out.stats.collect{it[1]}.ifEmpty([]))
        reports           = reports.mix(FASTQ_ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))

        // mix dna and rna in one channel
        bam_mapped = bam_mapped_dna.mix(bam_mapped_rna)

        // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
        // Except if and only if skipping markduplicates or saving mapped bams
        if (params.save_mapped || (params.skip_tools && params.skip_tools.split(',').contains('markduplicates'))) {

            // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta_with_fai)
            // Create CSV to restart from this step
            if (params.save_output_as_bam) CHANNEL_ALIGN_CREATE_CSV(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, params.outdir, params.save_output_as_bam)
            else CHANNEL_ALIGN_CREATE_CSV(BAM_TO_CRAM_MAPPING.out.cram.join(BAM_TO_CRAM_MAPPING.out.crai, failOnDuplicate: true, failOnMismatch: true), params.outdir, params.save_output_as_bam)

            // Gather used softwares versions
            versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
            versions = versions.mix(BAM_TO_CRAM_MAPPING.out.versions_samtools)
        }

        // Gather used softwares versions
        versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)
        versions = versions.mix(FASTQ_ALIGN.out.versions)

    }

    emit:
    bam_mapped_rna   = bam_mapped_rna  // second pass with RG tags
    bam_mapped_dna   = bam_mapped_dna  // second pass with RG tags
    bam_mapped       = bam_mapped      // for preprocessing
    cram_mapped      = cram_mapped     // for preprocessing
    reports          = reports
    versions         = versions
}
