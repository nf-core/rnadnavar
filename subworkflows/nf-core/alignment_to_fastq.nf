//
// BAM/CRAM to FASTQ conversion, paired end only
//
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_MAP_MAP      } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_UNMAP_UNMAP  } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_UNMAP_MAP    } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_MAP_UNMAP    } from '../../modules/nf-core/modules/samtools/view/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_UNMAP       } from '../../modules/nf-core/modules/samtools/merge/main'
include { SAMTOOLS_COLLATEFASTQ as COLLATE_FASTQ_UNMAP } from '../../modules/nf-core/modules/samtools/collatefastq/main'
include { SAMTOOLS_COLLATEFASTQ as COLLATE_FASTQ_MAP   } from '../../modules/nf-core/modules/samtools/collatefastq/main'
include { CAT_FASTQ                                    } from '../../modules/nf-core/modules/cat/fastq/main'

workflow ALIGNMENT_TO_FASTQ {
    take:
        input // channel: [meta, alignment (BAM or CRAM), index (optional)]
        fasta // optional: reference file if CRAM format and reference not in header

    main:
        ch_versions = Channel.empty()
        // extract mapped/unmapped reads
        // Differences between samtools view tasks come from modules.config

        // MAP - MAP [Extracts all reads that are paired and removes all unmapped] [-b -f1 -F12]
        SAMTOOLS_VIEW_MAP_MAP(input, fasta)
        // UNMAP - UNMAP [Extracts all unmapped reads/mates and removes all not primary alignment] [-b -f12 -F256]
        SAMTOOLS_VIEW_UNMAP_UNMAP(input, fasta)
        // UNMAP - MAP [Extracts all unmapped reads and removes all not primary alignment and unmapped mates] [-b -f4 -F264]
        SAMTOOLS_VIEW_UNMAP_MAP(input, fasta)
        // MAP - UNMAP [Extracts all unmapped mates and removes all not primary alignment and unmapped reads] [-b -f8 -F260]
        SAMTOOLS_VIEW_MAP_UNMAP(input, fasta)

        // Merge UNMAP
        all_unmapped_bam = SAMTOOLS_VIEW_UNMAP_UNMAP.out.bam
            .join(SAMTOOLS_VIEW_UNMAP_MAP.out.bam, remainder: true)
            .join(SAMTOOLS_VIEW_MAP_UNMAP.out.bam, remainder: true)
            .map{ meta, unmap_unmap, unmap_map, map_unmap ->
                [meta, [unmap_unmap, unmap_map, map_unmap]]
            }
        SAMTOOLS_MERGE_UNMAP(all_unmapped_bam, fasta)
        // Collate & convert unmapped
        COLLATE_FASTQ_UNMAP(SAMTOOLS_MERGE_UNMAP.out.bam)
        unmapped_reads = COLLATE_FASTQ_UNMAP.out.reads
            .map{ meta, reads_R1_R2, reads_other, reads_singleton ->
                [meta, reads_R1_R2]
            }

        // Collate & convert mapped
        COLLATE_FASTQ_MAP(SAMTOOLS_VIEW_MAP_MAP.out.bam)
        mapped_reads = COLLATE_FASTQ_MAP.out.reads
            .map{ meta, reads_R1_R2, reads_other, reads_singleton ->
                [meta, reads_R1_R2]
            }

        // join mapped & unmapped fastq
        reads_to_concat = mapped_reads.join(unmapped_reads)
            .map{ meta, mapped_reads, unmapped_reads ->
                [meta, [mapped_reads[0], mapped_reads[1], unmapped_reads[0], unmapped_reads[1]]]
            }

        // Concatenate Mapped_R1 with Unmapped_R1 and Mapped_R2 with Unmapped_R2
        CAT_FASTQ(reads_to_concat)

        // Gather versions of all tools used
        ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
        ch_versions = ch_versions.mix(COLLATE_FASTQ_MAP.out.versions)
        ch_versions = ch_versions.mix(COLLATE_FASTQ_UNMAP.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_MERGE_UNMAP.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MAP_MAP.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MAP_UNMAP.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_UNMAP_MAP.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_UNMAP_UNMAP.out.versions)

    emit:
        reads       = CAT_FASTQ.out.reads
        versions    = ch_versions
}