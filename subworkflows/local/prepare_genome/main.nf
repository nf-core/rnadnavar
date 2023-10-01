//
// PREPARE GENOME
//

// Initialize channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

include { BWA_INDEX as BWAMEM1_INDEX             } from '../../../modules/nf-core/bwa/index/main'
include { BWAMEM2_INDEX                          } from '../../../modules/nf-core/bwamem2/index/main'
include { DRAGMAP_HASHTABLE                      } from '../../../modules/nf-core/dragmap/hashtable/main'
include { GTF2BED                                } from '../../../modules/local/gtf2bed'                                       //addParams(options: params.genome_options)
include { GUNZIP as GUNZIP_GENE_BED              } from '../../../modules/nf-core/gunzip/main'                         //addParams(options: params.genome_options)
include { STAR_GENOMEGENERATE                    } from '../../../modules/nf-core/star/genomegenerate/main'            //addParams(options: params.star_index_options)
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                         } from '../../../modules/nf-core/samtools/faidx/main'
include { TABIX_TABIX as TABIX_DBSNP             } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_INDELS      } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_KNOWN_SNPS        } from '../../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as TABIX_PON               } from '../../../modules/nf-core/tabix/tabix/main'
include { HISAT2_EXTRACTSPLICESITES              } from '../../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                           } from '../../../modules/nf-core/hisat2/build/main'

workflow PREPARE_GENOME {
    take:
        dbsnp                   // channel: [optional]  dbsnp
        fasta                   // channel: [mandatory] fasta
        fasta_fai               // channel: [optional]  fasta_fai
        germline_resource       // channel: [optional]  germline_resource
        known_indels            // channel: [optional]  known_indels
        known_snps              // channel: [optional]  known_snps
        pon                     // channel: [optional]  pon


    main:

    fasta    = fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] }
    versions = Channel.empty()

    // If aligner is bwa-mem
    BWAMEM1_INDEX(fasta)     // If aligner is bwa-mem
    BWAMEM2_INDEX(fasta)     // If aligner is bwa-mem2
    DRAGMAP_HASHTABLE(fasta) // If aligner is dragmap


    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    SAMTOOLS_FAIDX(fasta, [['id':null], []])


    // the following are flattened and mapped in case the user supplies more than one value for the param
    // written for KNOWN_INDELS, but preemptively applied to the rest
    // [ file1, file2 ] becomes [ [ meta1, file1 ], [ meta2, file2 ] ]
    // outputs are collected to maintain a single channel for relevant TBI files
    TABIX_DBSNP(dbsnp.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map{ it -> [ [ id:it.baseName ], it ] })
    TABIX_KNOWN_SNPS(known_snps.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
    TABIX_KNOWN_INDELS(known_indels.flatten().map{ it -> [ [ id:it.baseName ], it ] } )
    TABIX_PON(pon.flatten().map{ it -> [ [ id:it.baseName ], it ] })

 //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    ch_gffread_version = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            GUNZIP_GTF (
                Channel.fromPath(params.gtf).map{ it -> [[id:it[0].baseName], it] }
            )
            ch_gtf = GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect()
            versions = versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.fromPath(params.gtf).collect()
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            GUNZIP_GFF (
                Channel.fromPath(params.gff).map{ it -> [[id:it[0].baseName], it] }
            )
            ch_gff = GUNZIP_GFF.out.gunzip.map{ meta, gff -> [gff] }.collect()
            versions = versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.fromPath(params.gff).collect()
        }

        GFFREAD (
            ch_gff
        )
        .gtf
        .set { ch_gtf }

        versions = versions.mix(GFFREAD.out.versions)
    }

     //
    // Uncompress exon BED annotation file or create from GTF if required
    //
    if (params.exon_bed) {
        if (params.exon_bed.endsWith('.gz')) {
            GUNZIP_GENE_BED (
                Channel.fromPath(params.exon_bed).map{ it -> [[id:it[0].baseName], it] }
            )
            ch_gene_bed = GUNZIP_GENE_BED.out.gunzip.map{ meta, bed -> [bed] }.collect()
            versions = versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_exon_bed = Channel.fromPath(params.exon_bed).collect()
        }
    } else {
        ch_exon_bed = GTF2BED ( ch_gtf ).bed.collect()
        versions = versions.mix(GTF2BED.out.versions)
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if (params.star_index) {
        if (params.star_index.endsWith('.tar.gz')) {
            UNTAR_STAR_INDEX (
                Channel.fromPath(params.star_index).map{ it -> [[id:it[0].baseName], it] }
            )
            ch_star_index = UNTAR_STAR_INDEX.out.untar.map{ meta, star_index -> [star_index] }.collect()
            versions   = versions.mix(UNTAR_STAR_INDEX.out.versions)
        } else {
            ch_star_index = Channel.fromPath(params.star_index).collect()
        }
    }
    else {
        STAR_GENOMEGENERATE ( fasta.map{meta, fasta -> fasta},ch_gtf )
        ch_star_index = STAR_GENOMEGENERATE.out.index
        versions      = versions.mix(STAR_GENOMEGENERATE.out.versions)
    }


    // HISAT2 not necessary if second pass skipped
    if ((params.skip_tools && !params.skip_tools.split(',').contains("second_run"))){
        if (params.splicesites) {
            ch_splicesites  = Channel.fromPath(params.splicesites).collect()
        } else{
            HISAT2_EXTRACTSPLICESITES ( ch_gtf.map{ it -> [ [ id:'null' ], it ]} )
            ch_splicesites  = HISAT2_EXTRACTSPLICESITES.out.txt
            versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        }

        if (params.hisat2_index) {
            ch_hisat2_index  = Channel.fromPath(params.hisat2_index).collect()
        } else{
			HISAT2_BUILD (
							fasta,
                            ch_gtf.map{ it -> [ [ id:'null' ], it ]},
                            ch_splicesites
                         )
            ch_hisat2_index = HISAT2_BUILD.out.index
            versions = versions.mix(HISAT2_BUILD.out.versions)
        }
    } else {
        ch_hisat2_index = Channel.empty()
        ch_splicesites  = Channel.empty()
    }


    // Gather versions of all tools used
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    versions = versions.mix(BWAMEM1_INDEX.out.versions)
    versions = versions.mix(BWAMEM2_INDEX.out.versions)
    versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    versions = versions.mix(TABIX_DBSNP.out.versions)
    versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    versions = versions.mix(TABIX_KNOWN_SNPS.out.versions)
    versions = versions.mix(TABIX_KNOWN_INDELS.out.versions)
    versions = versions.mix(TABIX_PON.out.versions)

    emit:
        bwa                   = BWAMEM1_INDEX.out.index.map{ meta, index -> [index] }.collect()       // path: bwa/*
        bwamem2               = BWAMEM2_INDEX.out.index.map{ meta, index -> [index] }.collect()       // path: bwamem2/*
        hashtable             = DRAGMAP_HASHTABLE.out.hashmap.map{ meta, index -> [index] }.collect() // path: dragmap/*
        dbsnp_tbi             = TABIX_DBSNP.out.tbi.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
        dict                  = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
        fasta_fai             = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
        germline_resource_tbi = TABIX_GERMLINE_RESOURCE.out.tbi.map{ meta, tbi -> [tbi] }.collect()   // path: germline_resource.vcf.gz.tbi
        known_snps_tbi        = TABIX_KNOWN_SNPS.out.tbi.map{ meta, tbi -> [tbi] }.collect()          // path: {known_indels*}.vcf.gz.tbi
        known_indels_tbi      = TABIX_KNOWN_INDELS.out.tbi.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        pon_tbi               = TABIX_PON.out.tbi.map{ meta, tbi -> [tbi] }.collect()                 // path: pon.vcf.gz.tbi
        star_index            = ch_star_index                                                         // path: star/index/
        gtf                   = ch_gtf                                                                // path: genome.gtf
        exon_bed              = ch_exon_bed                                                           // path: exon.bed
        hisat2_index          = ch_hisat2_index                                                       // path: hisat2/index/
        splicesites           = ch_splicesites                                                        // path: genome.splicesites.txt
        versions              = versions                                                              // channel: [ versions.yml ]
}
