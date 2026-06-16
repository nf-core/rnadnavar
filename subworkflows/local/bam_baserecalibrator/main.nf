//
// PREPARE RECALIBRATION
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'

workflow BAM_BASERECALIBRATOR {
    take:
    cram            // channel: [mandatory] [ meta, cram_markduplicates, crai ]
    dict            // channel: [mandatory] [ meta, dict ]
    fasta           // channel: [mandatory] [ meta, fasta ]
    fasta_fai       // channel: [mandatory] [ fasta_fai ]
    intervals       // channel: [mandatory] [ intervals, num_intervals ] (or [ [], 0 ] if no intervals)
    known_sites     // channel: [optional]  [ known_sites ]
    known_sites_tbi // channel: [optional]  [ known_sites_tbi ]

    main:
    versions = Channel.empty()
    known_sites_path = known_sites.map { entry -> entry instanceof List ? entry[-1] : entry }
    known_sites_tbi_path = known_sites_tbi.map { entry -> entry instanceof List ? entry[-1] : entry }

    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, c, idx, intervls, num_intervals -> [ meta + [ num_intervals:num_intervals ], c, idx, intervls ] }

    // RUN BASERECALIBRATOR
    GATK4_BASERECALIBRATOR(
                            cram_intervals,
                            fasta,
                            fasta_fai.map{ fai -> [ [id: "fai"], fai ] },
                            dict,
                            known_sites_path.map{ site -> [ [id: "sites"], site ] },
                            known_sites_tbi_path.map{ site_tbi -> [ [id: "sites_tbi"], site_tbi ] }
                            )

    // Figuring out if there is one or more table(s) from the same sample
    table_to_merge = GATK4_BASERECALIBRATOR.out.table.map{ meta, table -> [ groupKey(meta, meta.num_intervals), table ] }.groupTuple().branch{
        // Use meta.num_intervals to asses number of intervals
        single:   it[0].num_intervals <= 1
        multiple: it[0].num_intervals > 1
    }

    // Only when using intervals
    GATK4_GATHERBQSRREPORTS(table_to_merge.multiple)

    // Mix intervals and no_intervals channels together
    table_bqsr = GATK4_GATHERBQSRREPORTS.out.table.mix(table_to_merge.single.map{ meta, table -> [ meta, table[0] ] })
        // Remove no longer necessary field: num_intervals
        .map{ meta, table -> [ meta - meta.subMap('num_intervals'), table ] }

    // Gather versions of all tools used
    versions = versions.mix(GATK4_BASERECALIBRATOR.out.versions_gatk4)
    versions = versions.mix(GATK4_GATHERBQSRREPORTS.out.versions_gatk4)

    emit:
    table_bqsr // channel: [ meta, table ]

    versions   // channel: [ versions.yml ]
}
