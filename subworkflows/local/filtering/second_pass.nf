//
//
//
include { GATK_PREPROCESSING                                   } from '../gatk_preprocessing'
include { VARIANT_CALLING                                      } from '../variant_calling'
include { NORMALIZE                                            } from '../normalize_vcf_variants'
include { CONSENSUS                                            } from '../consensus'
include { GATK4_BEDTOINTERVALLIST                              } from '../../../modules/nf-core/modules/gatk4/bedtointervallist/main'
include { GATK4_FILTERSAMREADS                                 } from '../../../modules/nf-core/modules/gatk4/filtersamreads/main'
include { MAF2BED                                              } from '../../../modules/local/maf2bed'
//include { SAMTOOLS_CONVERT                                     } from '../../../modules/nf-core/modules/samtools/convert/main'
include { SAMTOOLS_INDEX                                       } from '../../../modules/nf-core/modules/samtools/index/main'
include { EXTRACT_READ_IDS                                     } from '../../../modules/local/extract_reads'
include { ALIGNMENT_TO_FASTQ as ALIGNMENT_TO_FASTQ_INPUT       } from '../../nf-core/alignment_to_fastq'
include { ALIGN_HISAT2                                         } from '../../nf-core/align_hisat2'
//include { SAMTOOLS_CONVERT as SAMTOOLS_BAMTOCRAM               } from '../../../modules/nf-core/modules/samtools/convert/main'
// Annotation
include { ANNOTATE                                             } from '../annotate'
include { VCF2MAF                                              } from '../../../modules/local/vcf2maf/vcf2maf/main'




workflow SECOND_PASS {
    take:
        ch_input_sample
        maf // for consensus
        bams // for realigment
        vcf_consensus_dna // to repeat rescue consensus
        vcfs_status_dna // to repeat rescue consensus
        fasta
        fasta_fai
        dict
        hisat2_index
        splicesites
        dbsnp
        dbsnp_tbi
        pon
        pon_tbi
        germline_resource
        germline_resource_tbi
        intervals
        intervals_for_preprocessing
        ch_interval_list_split
        intervals_bed_gz_tbi
        intervals_bed_combined

    main:
        ch_reports   = Channel.empty()
        ch_versions  = Channel.empty()
        // Extract allele coordinates from consensus

        MAF2BED(maf)
        bed = MAF2BED.out.bed
        bed.dump(tag:"[STEP8: RNA_FILTERING] bed1")
        // Match metadata for tumor and normal samples - this file comes from variant calling
        bed_normal = bed.map{meta, bed ->
                             [[
                id:         meta.normal_id + "_2pass",
                patient:    meta.patient,
                status:     0   // Had status from tumour because it comes from the variant calling consensus
                ],
            bed]}
        bed_tumor = bed.map{meta, bed ->
                             [[
                id:         meta.tumor_id + "_2pass" ,
                patient:    meta.patient,
                status:     meta.status
                ],
            bed]}

        bed = bed_normal.mix(bed_tumor).map{meta, bed -> [meta.id, meta, bed]}
        bed.dump(tag:"[STEP8: RNA_FILTERING] bed2")

        // Extract reads from BAMs
        SAMTOOLS_INDEX(bams) // index bams
        bams  = bams.join(SAMTOOLS_INDEX.out.bai)
        read_to_cross = bams.map{meta, bam, bai ->
                                                [ meta.id+ "_2pass",
                                                   [id:        meta.id + "_2pass",
                                                   patient:    meta.patient,
                                                   status:     meta.status,
                                                   read_group: meta.read_group],
                                                   bam,bai] }

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
        EXTRACT_READ_IDS(reads_to_filter, bed_to_filter)
        read_ids = EXTRACT_READ_IDS.out.read_ids
        read_ids.dump(tag:"[STEP8: RNA_FILTERING] read_ids")

        // Extract reads according to selected read ids
        bam_read = reads_to_filter.join(read_ids)
        GATK4_FILTERSAMREADS(bam_read, fasta) // bam -> filtered_bam
        GATK4_FILTERSAMREADS.out.bam.dump(tag:'[STEP8: RNA_FILTERING] filtered_bam')

        // Get FQs for re-alignment
        ALIGNMENT_TO_FASTQ_INPUT(GATK4_FILTERSAMREADS.out.bam, []) // bam -> fq
        ALIGNMENT_TO_FASTQ_INPUT.out.reads.dump(tag:'[STEP8: RNA_FILTERING] fastq_for_realignment')

        // HISAT2 re-alignment
        reads_to_realign = ALIGNMENT_TO_FASTQ_INPUT.out.reads.map{ meta, reads ->
                        // TODO: will this give a concurrent error too?
                        read_files = reads.sort{ a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                        [[
                            data_type:"fastq",
                            id:meta.id,
                            patient: meta.patient,
                            sample:meta.sample,
                            size:read_files.size(),
                            read_group:meta.read_group,
                            status:meta.status,
                            rg_id: meta.rg_id,
                            rg_pu: meta.rg_pu,
                            rg_sm: meta.rg_sm,
                            rg_lb: meta.rg_lb,
                            rg_pl: meta.rg_pl
                        ],
                        read_files[0]]
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
        ch_genome_bam_for_md = ch_genome_bam.map{meta, bam ->
                                        [[data_type:"bam",
                                        id:meta.id,
                                        patient: meta.patient,
                                        sample:meta.sample,
                                        status:meta.status],
                                        bam]
                                        }
                                   .groupTuple()
        ch_genome_bam.dump(tag:'[STEP8: RNA_FILTERING] ch_genome_bam')
        ch_genome_bam_index.dump(tag:'[STEP8: RNA_FILTERING] ch_genome_bam_index')

        // convert to CRAM -- no MD needs bams
//        SAMTOOLS_BAMTOCRAM(ch_genome_bam_bai, fasta, fasta_fai)
//        ch_genome_cram_crai = SAMTOOLS_BAMTOCRAM.out.alignment_index
        intervals.dump(tag:'[STEP8: RNA_FILTERING] intervals')
        intervals_for_preprocessing.dump(tag:'[STEP8: RNA_FILTERING] intervals_for_preprocessing')
        ch_interval_list_split.dump(tag:'[STEP8: RNA_FILTERING] ch_interval_list_split')
        ch_genome_bam_for_md.dump(tag:'[STEP8: RNA_FILTERING] ch_genome_bam_for_md')

        // STEP 2: GATK PREPROCESSING - See: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
    GATK_PREPROCESSING(
        'mapping',                   // Mandatory, step to start with
        ch_genome_bam_for_md,           // channel: [mandatory] [meta, [bam]]
        params.skip_tools,           // channel: [mandatory] skip_tools
        params.use_gatk_spark,       // channel: [mandatory] use_gatk_spark
        params.save_output_as_bam,   // channel: [mandatory] save_output_as_bam
        fasta,                       // channel: [mandatory] fasta
        fasta_fai ,                  // channel: [mandatory] fasta_fai
        dict,
        germline_resource,           // channel: [optional]  germline_resource
        germline_resource_tbi,       // channel: [optional]  germline_resource_tbi
        intervals,                   // channel: [mandatory] intervals/target regions
        intervals_for_preprocessing, // channel: [mandatory] intervals_for_preprocessing/wes
        ch_interval_list_split,
        ch_input_sample
    )

    ch_cram_variant_calling = GATK_PREPROCESSING.out.ch_cram_variant_calling
    ch_versions = ch_versions.mix(GATK_PREPROCESSING.out.versions)
    ch_reports = ch_reports.mix(GATK_PREPROCESSING.out.ch_reports)

//    intervals = intervals.map{meta, }
    ch_cram_variant_calling.dump(tag:"[STEP8 RNA_FILTERING] ch_cram_variant_calling")
    intervals_bed_gz_tbi.dump(tag:"[STEP8 RNA_FILTERING] intervals_bed_gz_tbi")
    pon.dump(tag:"[STEP8 RNA_FILTERING] pon")
// STEP 3: VARIANT CALLING
    VARIANT_CALLING(
        ch_cram_variant_calling,
        fasta,
        fasta_fai,
        dbsnp,
        dbsnp_tbi,
        dict,
        germline_resource,
        germline_resource_tbi,
        intervals,
        intervals_bed_gz_tbi,
        intervals_bed_combined,
        pon,
        pon_tbi
    )
    cram_vc_pair = VARIANT_CALLING.out.cram_vc_pair  // use same crams for force calling later
    vcf_to_normalize = VARIANT_CALLING.out.vcf
    ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)
    ch_reports = ch_reports.mix(VARIANT_CALLING.out.reports)


// STEP 4: NORMALIZE
   NORMALIZE (vcf_to_normalize,
              fasta )
   ch_versions = ch_versions.mix(NORMALIZE.out.versions)
   vcf_normalized = NORMALIZE.out.vcf

// STEP 5: CONSENSUS
    CONSENSUS (
                vcf_normalized,
                cram_vc_pair,  // from previous variant calling
                dict,
                fasta,
                fasta_fai,
                germline_resource,
                germline_resource_tbi,
                intervals,
                intervals_bed_gz_tbi,
                intervals_bed_combined,
                pon,
                pon_tbi,
                vcf_consensus_dna,
                vcfs_status_dna
                       )

// STEP 6: ANNOTATE
    ANNOTATE(
        ch_input_sample, // first pass
        CONSENSUS.out.vcf, // second pass TODO: make it optional
        fasta)

    ch_versions = ch_versions.mix(ANNOTATE.out.versions)
    ch_reports  = ch_reports.mix(ANNOTATE.out.reports)

    VCF2MAF(ANNOTATE.out.vcf_ann,
                fasta)
    maf_to_filter = VCF2MAF.out.maf

    emit:
        maf                             = maf_to_filter
        versions                        = ch_versions                                                         // channel: [ versions.yml ]
        reports                         = ch_reports

}