//
// PREPARE REFERENCE AND INTERVAL FILES FOR PIPELINE
//
include { PREPARE_GENOME                                       } from './prepare_genome'
include { PREPARE_INTERVALS                                    } from './prepare_intervals'
include { GATK4_BEDTOINTERVALLIST                              } from '../../modules/nf-core/modules/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS                              } from '../../modules/nf-core/modules/gatk4/intervallisttools/main'


workflow PREPARE_REFERENCE_AND_INTERVALS {

    main:
        ch_versions = Channel.empty()

        // Initialize file channels based on params, defined in the params.genomes[params.genome] scope
        dbsnp              = params.dbsnp              ? Channel.fromPath(params.dbsnp).collect()                    : Channel.value([])
        known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()               : Channel.value([])
        fasta              = params.fasta              ? Channel.fromPath(params.fasta).collect()                    : Channel.empty()
        fasta_fai          = params.fasta_fai          ? Channel.fromPath(params.fasta_fai).collect()                : Channel.empty()
        germline_resource  = params.germline_resource  ? Channel.fromPath(params.germline_resource).collect()        : Channel.value([]) //Mutec2 does not require a germline resource, so set to optional input
        known_indels       = params.known_indels       ? Channel.fromPath(params.known_indels).collect()             : Channel.value([])
        known_snps         = params.known_snps         ? Channel.fromPath(params.known_snps).collect()               : Channel.value([])
        pon                = params.pon                ? Channel.fromPath(params.pon).collect()                      : Channel.value([]) //PON is optional for Mutect2 (but highly recommended)
        whitelist          = params.whitelist          ? Channel.fromPath(params.whitelist).collect()                      : Channel.value([])

        // STEP 0.A: Build indices if needed
        PREPARE_GENOME(
            dbsnp,
            fasta,
            fasta_fai,
            germline_resource,
            known_indels,
            known_snps,
            pon)
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

        // Gather built indices or get them from the params
        bwa                    = params.fasta                   ? params.bwa                        ? Channel.fromPath(params.bwa).collect()                   : PREPARE_GENOME.out.bwa                   : []
        bwamem2                = params.fasta                   ? params.bwamem2                    ? Channel.fromPath(params.bwamem2).collect()               : PREPARE_GENOME.out.bwamem2               : []
        dragmap                = params.fasta                   ? params.dragmap                    ? Channel.fromPath(params.dragmap).collect()               : PREPARE_GENOME.out.hashtable             : []
        hisat2_index           = params.fasta                   ? params.hisat2_index               ? Channel.fromPath(params.hisat2_index).collect()          : PREPARE_GENOME.out.hisat2_index          : []
        splicesites            = params.fasta                   ? params.splicesites                ? Channel.fromPath(params.splicesites).collect()           : PREPARE_GENOME.out.splicesites           : []
        dict                   = params.fasta                   ? params.dict                       ? Channel.fromPath(params.dict).collect()                  : PREPARE_GENOME.out.dict                  : []
        fasta_fai              = params.fasta                   ? params.fasta_fai                  ? Channel.fromPath(params.fasta_fai).collect()             : PREPARE_GENOME.out.fasta_fai             : []
        dbsnp_tbi              = params.dbsnp                   ? params.dbsnp_tbi                  ? Channel.fromPath(params.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.value([])
        germline_resource_tbi  = params.germline_resource       ? params.germline_resource_tbi      ? Channel.fromPath(params.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : []
        known_indels_tbi       = params.known_indels            ? params.known_indels_tbi           ? Channel.fromPath(params.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.value([])
        known_snps_tbi         = params.known_snps              ? params.known_snps_tbi             ? Channel.fromPath(params.known_snps_tbi).collect()        : PREPARE_GENOME.out.known_snps_tbi        : Channel.value([])
        pon_tbi                = params.pon                     ? params.pon_tbi                    ? Channel.fromPath(params.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : []
        // known_sites is made by grouping both the dbsnp and the known snps/indels resources
        // Which can either or both be optional
        known_sites_indels     = dbsnp.concat(known_indels).collect()
        known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

        known_sites_snps     = dbsnp.concat(known_snps).collect()
        known_sites_snps_tbi = dbsnp_tbi.concat(known_snps_tbi).collect()

    // STEP 0.B: Build intervals if needed
        PREPARE_INTERVALS(fasta_fai)
        ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions)

        // Intervals for speed up preprocessing/variant calling by spread/gather
        intervals_bed_combined      = params.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined  // [interval.bed] all intervals in one file
        intervals_for_preprocessing = params.wes          ? intervals_bed_combined : []      // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
        intervals                   = PREPARE_INTERVALS.out.intervals_bed        // [interval, num_intervals] multiple interval.bed files, divided by useful intervals for scatter/gather
        intervals_bed_gz_tbi        = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [interval_bed, tbi, num_intervals] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather


    // STEP 0.C: Prepare the interval list from the GTF file using GATK4 BedToIntervalList
        ch_genome_bed = Channel.from([id:'genome.bed']).combine(PREPARE_GENOME.out.exon_bed)
        ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)
        ch_interval_list = Channel.empty()
        GATK4_BEDTOINTERVALLIST(
            ch_genome_bed,
            dict
            )
        ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
        ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

    // STEP 0.D: Scatter one interval-list into many interval-files using GATK4 IntervalListTools
        ch_interval_list_split = Channel.empty()
        if (!params.skip_intervallisttools) {
            GATK4_INTERVALLISTTOOLS(
                ch_interval_list
            )
            ch_interval_list_split = GATK4_INTERVALLISTTOOLS.out.interval_list.map{ meta, bed -> [bed] }.flatten()
        }
        else {
            ch_interval_list_split = ch_interval_list
        }

    emit:
        fasta = fasta
        fasta_fai = fasta_fai
        dict = dict
        bwa = bwa
        germline_resource = germline_resource
        germline_resource_tbi = germline_resource_tbi
        bwamem2 = bwamem2
        dragmap                     = dragmap
        star_index                  = PREPARE_GENOME.out.star_index
        gtf                         = PREPARE_GENOME.out.gtf
        ch_interval_list            = ch_interval_list
        ch_interval_list_split      = ch_interval_list_split
        intervals                   = intervals
        intervals_bed_gz_tbi        = intervals_bed_gz_tbi
        intervals_for_preprocessing = intervals_for_preprocessing
        intervals_bed_combined      = intervals_bed_combined
        dbsnp                       = dbsnp
        dbsnp_tbi                   = dbsnp_tbi
        pon                         = pon
        pon_tbi                     = pon_tbi
        germline_resource           = germline_resource
        germline_resource_tbi       = germline_resource_tbi
        hisat2_index                = hisat2_index
        splicesites                 = splicesites
        versions                    = ch_versions                                            // channel: [ versions.yml ]

}