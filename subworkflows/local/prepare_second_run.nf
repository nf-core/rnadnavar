//
//
//

include { BAM_MERGE_INDEX_SAMTOOLS as MERGE_ALIGN     } from './bam_merge_index_samtools/main'
include { MAF2BED                                     } from '../../modules/local/maf2bed'
include { EXTRACT_READ_IDS                            } from '../../modules/local/extract_reads'
include { ALIGNMENT_TO_FASTQ as ALIGN2FQ              } from '../nf-core/alignment_to_fastq'
include { ALIGN_HISAT2                                } from '../nf-core/align_hisat2'
include { GATK4_FILTERSAMREADS                        } from '../../modules/nf-core/modules/gatk4/filtersamreads/main'


workflow PREPARE_SECOND_RUN {
    take:
        ch_input_sample
        maf // for consensus
        bwa_bams // for realigment
        star_bams // for realigment
        fasta
        fasta_fai
        dict
        hisat2_index
        splicesites

    main:
        ch_reports   = Channel.empty()
        ch_versions  = Channel.empty()

        if (params.step in ['mapping', 'markduplicates', 'splitncigar', 'prepare_recalibration', 'recalibrate', 'variant_calling', 'normalize', 'consensus', 'second_pass'] ) {

            // RNA specific filtering (2nd PASS) - this is fast BUT it increases the length of the pipeline considerably
            star_bams.dump(tag:"star_bams")
            bwa_bams.dump(tag:"bwa_bams")
            // 1 We take the previous aligned reads with star for tumor RNA and DNA normal
            bwa_bams.branch{
                            normal: it[0].status == 0
                            tumor: it[0].status == 1
                        }.set{previous_dna_alignment}
            // we only need normals - dna tumour will NOT be realigned
            previous_normal_alignment = previous_dna_alignment.normal.groupTuple()
            // 2. Group them and merge if applicable
            previous_alignment = star_bams.map{meta, bam, bai -> [meta, [bam], bai]}
                                        .mix(previous_normal_alignment)
            previous_alignment.dump(tag: 'previous_alignment')
            MERGE_ALIGN(previous_alignment.map{meta, bam, bai -> [meta, bam]})
            ch_versions = ch_versions.mix(MERGE_ALIGN.out.versions)
            // TODO make indexing optional if bai already there
            MERGE_ALIGN.out.bam_bai.dump(tag:'MERGE_ALIGN.out.bam_bai')
            previous_bams = MERGE_ALIGN.out.bam_bai
                                            .map{meta, bam, bai -> [
                                                                [
                                                                    id:meta.sample,
                                                                    data_type:"bam",
                                                                    patient:meta.patient,
                                                                    sample:meta.sample,
                                                                    read_group:meta.read_group,
                                                                    status:meta.status
                                                                ],
                                                                bam, bai
                                                              ]
                                                 }
            previous_bams.dump(tag:"[STEP7: FILTERING] bams for realignment")

    // STEP A: Extract allele coordinates from consensus
            MAF2BED(maf)
            bed = MAF2BED.out.bed
            bed.branch{
                       dna: it[0].status == 1  // Bed from DNA tumour
                       rna: it[0].status == 2   // Bed from RNA tumour
            }.set{bed_status}
            bed.dump(tag:"[STEP8: RNA_FILTERING] bed1")
            // Match metadata for tumor and normal samples - this file comes from variant calling
            bed_normal = bed_status.rna.map{meta, bed ->
                    def (tumor_id, normal_id) = meta.rna_id.tokenize( '_vs_' )
                                 [[
                    id:         normal_id + "_2pass",
                    patient:    meta.patient,
                    status:     0   // Had status from tumour because it comes from the variant calling consensus
                    ],
                bed]}
            bed_tumor = bed_status.rna.map{meta, bed ->
                    def (tumor_id, normal_id) = meta.rna_id.tokenize( '_vs_' )
                                 [[
                    id:         tumor_id + "_2pass" ,
                    patient:    meta.patient,
                    status:     meta.status
                    ],
                bed]}
            bed_normal.dump(tag:"bed_normal")
            bed_tumor.dump(tag:"bed_tumor")
            bed = bed_normal.mix(bed_tumor).map{meta, bed -> [meta.id, meta, bed]}
            bed.dump(tag:"[STEP8: RNA_FILTERING] bed2")

    // STEP B: Extract reads from BAMs
            read_to_cross = previous_bams.map{meta, bam, bai ->
                                                    [ meta.id+ "_2pass",
                                                       [id:        meta.id + "_2pass",
                                                       patient:    meta.patient,
                                                       status:     meta.status,
                                                       read_group: meta.read_group],
                                                       bam,bai] }
            read_to_cross.dump(tag:'read_to_cross')

            crossed_reads_bed = read_to_cross.cross(bed)
            // We cross the values to make sure we match per id the bed file with the bams
            bed_to_filter = crossed_reads_bed.map{ reads, bed ->
                                                     def rg_id = ( reads[1].read_group =~ /ID:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_pu = ( reads[1].read_group =~ /PU:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_sm = ( reads[1].read_group =~ /SM:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_lb = ( reads[1].read_group =~ /LB:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_pl = ( reads[1].read_group =~ /PL:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def meta = [:]
                                                      meta.id         = bed[1].id
                                                      meta.patient    = bed[1].patient
                                                      meta.status     = bed[1].status
                                                      meta.sample     = rg_sm + "_2pass"
                                                      meta.read_group = reads[1].read_group
                                                      meta.rg_id = rg_id
                                                      meta.rg_pu = rg_pu
                                                      meta.rg_sm = rg_sm + "_2pass"
                                                      meta.rg_lb = rg_lb
                                                      meta.rg_pl = rg_pl
                                                      [meta, bed[2]]}
            reads_to_filter = crossed_reads_bed.map{ reads, bed ->
                                                     def rg_id = ( reads[1].read_group =~ /ID:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_pu = ( reads[1].read_group =~ /PU:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_sm = ( reads[1].read_group =~ /SM:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_lb = ( reads[1].read_group =~ /LB:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def rg_pl = ( reads[1].read_group =~ /PL:(\S+?)(\\t|\s|\"|$)/ )[0][1]
                                                     def meta = [:]
                                                      meta.id         = reads[1].id
                                                      meta.patient    = reads[1].patient
                                                      meta.status     = reads[1].status
                                                      meta.sample     = rg_sm + "_2pass"
                                                      meta.read_group = reads[1].read_group
                                                      meta.rg_id = rg_id
                                                      meta.rg_pu = rg_pu
                                                      meta.rg_sm = rg_sm + "_2pass"
                                                      meta.rg_lb = rg_lb
                                                      meta.rg_pl = rg_pl
                                                      [meta, reads[2], reads[3]]}
            reads_to_filter.dump(tag:"[STEP8: RNA_FILTERING] reads_to_filter")
            bed_to_filter.dump(tag:"[STEP8: RNA_FILTERING] bed_to_filter")
            EXTRACT_READ_IDS(reads_to_filter, bed_to_filter)
            read_ids = EXTRACT_READ_IDS.out.read_ids
            read_ids.dump(tag:"[STEP8: RNA_FILTERING] read_ids")

    //STEP C: Extract reads according to selected read ids
            bam_read = reads_to_filter.join(read_ids)
            GATK4_FILTERSAMREADS(bam_read, fasta) // bam -> filtered_bam
            GATK4_FILTERSAMREADS.out.bam.dump(tag:'[STEP8: RNA_FILTERING] filtered_bam')

    //STEP D: Get FQs for re-alignment
            ALIGN2FQ(GATK4_FILTERSAMREADS.out.bam, []) // bam -> fq
            ALIGN2FQ.out.reads.dump(tag:'[STEP8: RNA_FILTERING] fastq_for_realignment')

    //STEP E: HISAT2 re-alignment
            reads_to_realign = ALIGN2FQ.out.reads.map{ meta, reads ->
                            // This can throw a concurrent error when submitting a lot of samples
//                            read_files = reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                            [[data_type:"fastq"] + meta, reads]
                            }
            reads_to_realign.dump(tag:'[STEP8: RNA_FILTERING] reads_to_realign')
            // HISAT2 realignment
            ALIGN_HISAT2 (
                reads_to_realign,
                hisat2_index,
                splicesites
            )
            ch_genome_bam        = ALIGN_HISAT2.out.bam
            ch_genome_bam_index  = ALIGN_HISAT2.out.bai
            ch_samtools_stats    = ALIGN_HISAT2.out.stats
            ch_samtools_flagstat = ALIGN_HISAT2.out.flagstat
            ch_samtools_idxstats = ALIGN_HISAT2.out.idxstats
            ch_hisat2_multiqc    = ALIGN_HISAT2.out.summary
            ch_genome_bam.join(ch_genome_bam_index).dump(tag:"HERE")
            ch_genome_bam.join(ch_genome_bam_index).map{meta, bam, bai ->
                                            [[data_type:"bam",
                                            id:meta.id,
                                            patient: meta.patient,
                                            sample:meta.sample,
                                            status:meta.status],
                                            bam, bai]
                                            }.dump(tag:"HERE2")
            ch_genome_bam_for_md = ch_genome_bam.join(ch_genome_bam_index)
                                                .map{meta, bam, bai ->
                                                        [[data_type:"bam",
                                                        id:meta.id,
                                                        patient: meta.patient,
                                                        sample:meta.sample,
                                                        status:meta.status],
                                                        bam, bai]
                                                    }
            ch_genome_bam.dump(tag:'[STEP8: RNA_FILTERING] ch_genome_bam')
            ch_genome_bam_index.dump(tag:'[STEP8: RNA_FILTERING] ch_genome_bam_index')
//            ch_genome_bam_for_md = ch_genome_bam.mix(ch_genome_bam_index)

            ch_genome_bam_for_md.dump(tag:'[STEP8: RNA_FILTERING] ch_genome_bam_for_md')
    }

    emit:
        ch_bam_mapped    = ch_genome_bam_for_md
        versions         = ch_versions                                                         // channel: [ versions.yml ]
        reports          = ch_reports








}