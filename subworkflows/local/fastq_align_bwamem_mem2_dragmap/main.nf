//
// MAPPING
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BWAMEM2_MEM            } from '../../../modules/nf-core/bwamem2/mem/main'
include { BWA_MEM as BWAMEM1_MEM } from '../../../modules/nf-core/bwa/mem/main'
include { DRAGMAP_ALIGN          } from '../../../modules/nf-core/dragmap/align/main'

workflow FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP {
    take:
    reads // channel: [mandatory] meta, reads
    index // channel: [mandatory] meta, index
    fasta // channel: [mandatory] meta, fasta
    sort  // boolean: [mandatory] true -> sort, false -> don't sort

    main:

    versions = Channel.empty()
    reports = Channel.empty()

    // Only one of the following should be run
    BWAMEM1_MEM(reads, index.map{ it -> [ [ id:'index' ], it ] }, fasta, sort) // If aligner is bwa-mem
    BWAMEM2_MEM(reads, index.map{ it -> [ [ id:'index' ], it ] }, fasta, sort) // If aligner is bwa-mem2
    DRAGMAP_ALIGN(reads, index.map{ it -> [ [ id:'index' ], it ] }, fasta, sort) // If aligner is dragmap

    // Get the bam files from the aligner
    // Only one aligner is run
    bam = Channel.empty()
    bam = bam.mix(BWAMEM1_MEM.out.bam)
    bam = bam.mix(BWAMEM2_MEM.out.bam)
    bam = bam.mix(DRAGMAP_ALIGN.out.bam)

    // Gather reports of all tools used
    reports = reports.mix(DRAGMAP_ALIGN.out.log)

    // Gather versions of all tools used
    versions = versions.mix(BWAMEM1_MEM.out.versions)
    versions = versions.mix(BWAMEM2_MEM.out.versions)
    versions = versions.mix(DRAGMAP_ALIGN.out.versions)

    emit:
    bam      // channel: [ [meta], bam ]
    reports
    versions // channel: [ versions.yml ]
}
