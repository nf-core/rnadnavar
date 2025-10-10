//
// RECALIBRATE
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_APPLYBQSR           } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { CRAM_MERGE_INDEX_SAMTOOLS } from '../cram_merge_index_samtools/main'

workflow BAM_APPLYBQSR {
    take:
    cram          // channel: [mandatory] [ meta, cram, crai, recal ]
    dict          // channel: [mandatory] [ meta, dict ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] [ fasta_fai ]
    intervals     // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, c, crai, recal, intervls, num_intervals -> [ meta + [ num_intervals:num_intervals ], c, crai, recal, intervls ] }
    cram_intervals.dump(tag:"cram_intervalsGATK4_APPLYBQSR")
    // RUN APPLYBQSR
    GATK4_APPLYBQSR(
                    cram_intervals,
                    fasta,
                    fasta_fai,
                    dict.map{ _meta, it -> [ it ] }
                    )
    GATK4_APPLYBQSR.out.cram.dump(tag:"GATK4_APPLYBQSR.out.cram")
    // Gather the recalibrated cram files
    cram_to_merge = GATK4_APPLYBQSR.out.cram.map{ meta, c -> [ groupKey(meta, meta.num_intervals), c ] }.groupTuple()
    cram_to_merge.dump(tag:"cram_to_merge1")
    // Merge and index the recalibrated cram files
    CRAM_MERGE_INDEX_SAMTOOLS(cram_to_merge, fasta.map{_meta, fa -> fa}, fasta_fai)

    cram_recal = CRAM_MERGE_INDEX_SAMTOOLS.out.cram_crai
        // Remove no longer necessary field: num_intervals
        .map{ meta, c, idx -> [ meta - meta.subMap('num_intervals'), c, idx ] }

    // Gather versions of all tools used
    versions = versions.mix(GATK4_APPLYBQSR.out.versions)
    versions = versions.mix(CRAM_MERGE_INDEX_SAMTOOLS.out.versions)

    emit:
    cram = cram_recal // channel: [ meta, cram, crai ]

    versions          // channel: [ versions.yml ]
}
