//
// PREPARE INTERVALS
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BUILD_INTERVALS                                        } from '../../../modules/local/build_intervals/main'
include { CREATE_INTERVALS_BED                                   } from '../../../modules/local/create_intervals_bed/main'
include { GATK4_INTERVALLISTTOBED                                } from '../../../modules/nf-core/gatk4/intervallisttobed/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_SPLIT    } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_INTERVAL_COMBINED } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow PREPARE_INTERVALS {
    take:
    fasta_fai    // mandatory [ fasta_fai ]
    intervals    // [ params.intervals ]
    no_intervals // [ params.no_intervals ]

    main:
    versions = Channel.empty()


    if (no_intervals) {
        intervals_bed                          = Channel.of([[], 0 ])
        intervals_bed_gz_tbi                   = Channel.of([[[],[]], 0 ])
        intervals_combined                     = Channel.of([[id:"no_intervals"], 0 ])
        intervals_bed_gz_tbi_combined          = Channel.of([[], []])
        intervals_bed_gz_tbi_and_num_intervals = Channel.of([[],[], 0 ])
    } else if (params.step != 'annotate') {
        // If no interval/target file is provided, then generated intervals from FASTA file
        if (!intervals) {
            BUILD_INTERVALS(fasta_fai.map{it -> [ [ id:it.baseName ], it ] })

            intervals_combined = BUILD_INTERVALS.out.bed

            CREATE_INTERVALS_BED(intervals_combined.map{ meta, path -> path }).bed

            intervals_bed = CREATE_INTERVALS_BED.out.bed

            versions = versions.mix(BUILD_INTERVALS.out.versions)
            versions = versions.mix(CREATE_INTERVALS_BED.out.versions)
        } else {
            intervals_combined = Channel.fromPath(file(intervals)).map{it -> [ [ id:it.baseName ], it ] }
            intervals_bed = CREATE_INTERVALS_BED(file(intervals)).bed

            versions = versions.mix(CREATE_INTERVALS_BED.out.versions)

            // If interval file is not provided as .bed, but e.g. as .interval_list then convert to BED format
            if (intervals.endsWith(".interval_list")) {
                GATK4_INTERVALLISTTOBED(intervals_combined)
                intervals_combined = GATK4_INTERVALLISTTOBED.out.bed
                versions = versions.mix(GATK4_INTERVALLISTTOBED.out.versions)
            }
        }

        // Now for the intervals.bed the following operations are done:
        // 1. Intervals file is split up into multiple bed files for scatter/gather
        // 2. Each bed file is indexed

        // 1. Intervals file is split up into multiple bed files for scatter/gather & grouping together small intervals
        intervals_bed = intervals_bed.flatten()
            .map{ intervalFile ->
                def duration = 0.0
                for (line in intervalFile.readLines()) {
                    final fields = line.split('\t')
                    if (fields.size() >= 5) duration += fields[4].toFloat()
                    else {
                        start = fields[1].toInteger()
                        end = fields[2].toInteger()
                        duration += (end - start) / params.nucleotides_per_second
                    }
                }
                [ duration, intervalFile ]
            }.toSortedList({ a, b -> b[0] <=> a[0] })
            .flatten().collate(2).map{ duration, intervalFile -> intervalFile }.collect()
            // Adding number of intervals as elements
            .map{ it -> [ it, it.size() ] }
            .transpose()

        // 2. Create bed.gz and bed.gz.tbi for each interval file. They are split by region (see above)
        TABIX_BGZIPTABIX_INTERVAL_SPLIT(intervals_bed.map{ file, num_intervals -> [ [ id:file.baseName], file ] })

        intervals_bed_gz_tbi = TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.gz_tbi.map{ meta, bed, tbi -> [ bed, tbi ] }.toList()
            // Adding number of intervals as elements
            .map{ it -> [ it, it.size() ] }
            .transpose()

        versions = versions.mix(TABIX_BGZIPTABIX_INTERVAL_SPLIT.out.versions)

        TABIX_BGZIPTABIX_INTERVAL_COMBINED(intervals_combined)
        versions = versions.mix(TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.versions)


        intervals_bed_gz_tbi_combined = TABIX_BGZIPTABIX_INTERVAL_COMBINED.out.gz_tbi.map{meta, gz, tbi -> [gz, tbi] }.collect()
        intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map{ intervals, num_intervals ->
        if ( num_intervals < 1 ) [ [], [], num_intervals ]
        else [ intervals[0], intervals[1], num_intervals ]
        }
    }

    // Intervals for speed up preprocessing/variant calling by spread/gather
    intervals_bed_combined = no_intervals ?
        Channel.value([]) :
        intervals_combined.map{meta, bed -> bed }.collect()
    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = params.wes ?
        intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect() :
        Channel.value([ [ id:'null' ], [] ])

    intervals_and_num_intervals   = intervals_bed.map{ interval, num_intervals ->
            if ( num_intervals < 1 ) [ [], num_intervals ]
            else [ interval, num_intervals ]
        }

    emit:
    // Intervals split for parallel execution
    intervals_bed                 // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi          // [interval_bed, tbi, num_intervals] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather
    intervals_for_preprocessing
    // All intervals in one file
    intervals_bed_combined        // [ intervals.bed ]
    intervals_bed_gz_tbi_combined // [ intervals.bed.gz, intervals.bed.gz.tbi]
    intervals_and_num_intervals
    intervals_bed_gz_tbi_and_num_intervals

    versions               // [ versions.yml ]
}
