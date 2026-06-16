//
// PREPARE REFERENCE AND INTERVAL FILES FOR PIPELINE
//
include { PREPARE_GENOME                                       } from './../prepare_genome/main'
include { PREPARE_INTERVALS                                    } from './../prepare_intervals/main'
include { GATK4_BEDTOINTERVALLIST                              } from '../../../modules/nf-core/gatk4/bedtointervallist/main'

workflow PREPARE_REFERENCE_AND_INTERVALS {
    take:
    dbsnp_input
    known_snps_input
    fasta_input
    germline_resource_input
    known_indels_input
    pon_input
    whitelist_input
    bwa_input
    bwamem2_input
    dragmap_input
    hisat2_index_input
    splicesites_input
    dict_input
    fasta_fai_input
    dbsnp_tbi_input
    germline_resource_tbi_input
    known_indels_tbi_input
    known_snps_tbi_input
    pon_tbi_input
    intervals_input
    no_intervals_input
    step
    nucleotides_per_second
    wes

    main:
    versions = Channel.empty()

    dbsnp              = dbsnp_input              ? Channel.fromPath(dbsnp_input).collect()                    : Channel.value([])
    known_snps         = known_snps_input         ? Channel.fromPath(known_snps_input).collect()               : Channel.value([])
    fasta              = fasta_input              ? Channel.fromPath(fasta_input).collect()                    : Channel.empty()
    germline_resource  = germline_resource_input  ? Channel.fromPath(germline_resource_input).collect()        : Channel.value([]) //Mutec2 does not require a germline resource, so set to optional input
    known_indels       = known_indels_input       ? Channel.fromPath(known_indels_input).collect()             : Channel.value([])
    pon                = pon_input                ? Channel.fromPath(pon_input).collect()                      : Channel.value([]) //PON is optional for Mutect2 (but highly recommended)
    whitelist          = whitelist_input          ? Channel.fromPath(whitelist_input).collect()                : Channel.value([]) // whitelist optional for filtering

    // Build indexes if needed
    PREPARE_GENOME(
        dbsnp,
        fasta,
        germline_resource,
        known_indels,
        known_snps,
        pon)
    versions = versions.mix(PREPARE_GENOME.out.versions)

    // Gather built indices or get them from the params
    fasta                  = PREPARE_GENOME.out.fasta   // unbgzipped if .gz
    bwa                    = fasta_input                   ? bwa_input                        ? Channel.fromPath(bwa_input).collect()                           : PREPARE_GENOME.out.bwa                   : []
    bwamem2                = fasta_input                   ? bwamem2_input                    ? Channel.fromPath(bwamem2_input).collect()                       : PREPARE_GENOME.out.bwamem2               : []
    dragmap                = fasta_input                   ? dragmap_input                    ? Channel.fromPath(dragmap_input).collect()                       : PREPARE_GENOME.out.hashtable             : []
    hisat2_index           = fasta_input                   ? hisat2_index_input               ? Channel.fromPath(hisat2_index_input).map{ it -> [ [id:'ht_idx'], it ] }.collect()                          : PREPARE_GENOME.out.hisat2_index : []
    splicesites            = fasta_input                   ? splicesites_input                ? Channel.fromPath(splicesites_input).collect()                   : PREPARE_GENOME.out.splicesites           : []
    dict                   = dict_input                    ? Channel.fromPath(dict_input).map{ it -> [ [id:'dict'], it ] }.collect()                             : PREPARE_GENOME.out.dict
    fasta_fai              = fasta_input                   ? fasta_fai_input                  ? Channel.fromPath(fasta_fai_input).collect()                     : PREPARE_GENOME.out.fasta_fai             : []
    // If no external index was supplied, consume the rebuilt resource emitted by PREPARE_GENOME
    // so the resource path and generated TBI keep matching basenames downstream.
    dbsnp                  = dbsnp_input                   ? dbsnp_tbi_input                  ? Channel.fromPath(dbsnp_input).collect()                         : PREPARE_GENOME.out.dbsnp                  : Channel.value([])
    germline_resource      = germline_resource_input       ? germline_resource_tbi_input      ? Channel.fromPath(germline_resource_input).collect()             : PREPARE_GENOME.out.germline_resource      : []
    known_indels           = known_indels_input            ? known_indels_tbi_input           ? Channel.fromPath(known_indels_input).collect()                  : PREPARE_GENOME.out.known_indels           : Channel.value([])
    known_snps             = known_snps_input              ? known_snps_tbi_input             ? Channel.fromPath(known_snps_input).collect()                    : PREPARE_GENOME.out.known_snps             : Channel.value([])
    pon                    = pon_input                     ? pon_tbi_input                    ? Channel.fromPath(pon_input).collect()                           : PREPARE_GENOME.out.pon                    : []
    dbsnp_tbi              = dbsnp_input                   ? dbsnp_tbi_input                  ? Channel.fromPath(dbsnp_tbi_input).collect()                     : PREPARE_GENOME.out.dbsnp_tbi             : Channel.value([])
    germline_resource_tbi  = germline_resource_input       ? germline_resource_tbi_input      ? Channel.fromPath(germline_resource_tbi_input).collect()         : PREPARE_GENOME.out.germline_resource_tbi : []
    known_indels_tbi       = known_indels_input            ? known_indels_tbi_input           ? Channel.fromPath(known_indels_tbi_input).collect()              : PREPARE_GENOME.out.known_indels_tbi      : Channel.value([])
    known_snps_tbi         = known_snps_input              ? known_snps_tbi_input             ? Channel.fromPath(known_snps_tbi_input).collect()                : PREPARE_GENOME.out.known_snps_tbi        : Channel.value([])
    pon_tbi                = pon_input                     ? pon_tbi_input                    ? Channel.fromPath(pon_tbi_input).collect()                       : PREPARE_GENOME.out.pon_tbi               : []
    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    known_sites_snps       = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi   = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(
        fasta_fai,
        intervals_input,
        no_intervals_input,
        step,
        nucleotides_per_second,
        wes
    )
    versions = versions.mix(PREPARE_INTERVALS.out.versions)


    emit:
    fasta
    fasta_fai
    dict
    bwa
    germline_resource
    germline_resource_tbi
    bwamem2
    dragmap
    star_index                    = PREPARE_GENOME.out.star_index
    gtf                           = PREPARE_GENOME.out.gtf
    intervals                     = PREPARE_INTERVALS.out.intervals_bed
    intervals_bed_gz_tbi          = PREPARE_INTERVALS.out.intervals_bed_gz_tbi
    intervals_for_preprocessing   = PREPARE_INTERVALS.out.intervals_for_preprocessing
    intervals_bed_combined        = PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined = PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined
    dbsnp
    dbsnp_tbi
    pon
    pon_tbi
    hisat2_index
    splicesites
    known_sites_indels
    known_sites_indels_tbi
    known_sites_snps
    known_sites_snps_tbi
    versions                                                                  // channel: [ versions.yml ]
}
