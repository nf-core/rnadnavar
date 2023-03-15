//
// Subworkflow: Run GATK4 SplitNCigarReads without intervals, merge and index BAM file.
//
include { GATK4_SPLITNCIGARREADS } from '../../modules/nf-core/modules/gatk4/splitncigarreads/main'
include { SAMTOOLS_INDEX         } from '../../modules/nf-core/modules/samtools/index/main'

workflow SPLITNCIGAR {
    take:
        bam             // channel: [ val(meta), [ bam ], [bai] ]
        fasta           // channel: [ fasta ]
        fasta_fai       // channel: [ fai ]
        fasta_dict      // channel: [ dict ]
        intervals       // channel: [ interval_list]

    main:

        ch_versions       = Channel.empty()
        bam = bam.map{meta, bam, bai -> [meta, bam, bai, []]}
        GATK4_SPLITNCIGARREADS (
            bam,
            fasta,
            fasta_fai,
            fasta_dict
        )
        bam_splitncigar = GATK4_SPLITNCIGARREADS.out.bam
        ch_versions = ch_versions.mix(GATK4_SPLITNCIGARREADS.out.versions.first())

        SAMTOOLS_INDEX (
            bam_splitncigar
        )
        splitncigar_bam_bai = bam_splitncigar
            .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
            .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
            .map{meta, bam, bai, csi ->
                if (bai) [meta, bam, bai]
                else [meta, bam, csi]
            }
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
        bam_bai     = splitncigar_bam_bai
        versions    = ch_versions
}
