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
include { GUNZIP as GUNZIP_GENE_BED              } from '../../../modules/nf-core/gunzip/main'                         //addParams(options: params.genome_options)
include { GUNZIP as GUNZIP_GTF                   } from '../../../modules/nf-core/gunzip/main'                         //addParams(options: params.genome_options)
include { GUNZIP as GUNZIP_GFF                   } from '../../../modules/nf-core/gunzip/main'
include { GFFREAD                                } from '../../../modules/nf-core/gffread/main'
include { UNTAR as UNTAR_STAR_INDEX              } from '../../../modules/nf-core/untar/main'
include { TABIX_BGZIP as UNBGZIP_FASTA           } from '../../../modules/nf-core/tabix/bgzip/main'
include { HTSLIB_BGZIPTABIX as TABIX_DBSNP       } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { HTSLIB_BGZIPTABIX as TABIX_GERMLINE_RESOURCE } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { HTSLIB_BGZIPTABIX as TABIX_KNOWN_INDELS } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { HTSLIB_BGZIPTABIX as TABIX_KNOWN_SNPS  } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { HTSLIB_BGZIPTABIX as TABIX_PON         } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { STAR_GENOMEGENERATE                    } from '../../../modules/nf-core/star/genomegenerate/main'            //addParams(options: params.star_index_options)
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../../../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                         } from '../../../modules/nf-core/samtools/faidx/main'
include { HISAT2_EXTRACTSPLICESITES              } from '../../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                           } from '../../../modules/nf-core/hisat2/build/main'

workflow PREPARE_GENOME {
    take:
        dbsnp                   // channel: [optional]  dbsnp
        fasta                   // channel: [mandatory] fasta
        germline_resource       // channel: [optional]  germline_resource
        known_indels            // channel: [optional]  known_indels
        known_snps              // channel: [optional]  known_snps
        pon                     // channel: [optional]  pon

    main:


    versions = Channel.empty()

    // UNBGZIP genome if applicable
    if (params.fasta.endsWith('.gz')){  // bgzip
        UNBGZIP_FASTA ( Channel.fromPath(params.fasta).collect().map{ fa -> [ [ id:fa.baseName[0] - ~/\.fa(sta)?$/ ], fa ] } )
        fasta    = UNBGZIP_FASTA.out.output
        versions = versions.mix(UNBGZIP_FASTA.out.versions_tabix)
    } else{
        fasta    = fasta.map{ fa -> [ [ id:fa.baseName ], fa ] }
    }

    // If aligner is bwa-mem
    if (params.dna){
        BWAMEM1_INDEX(fasta)     // If aligner is bwa-mem
        BWAMEM2_INDEX(fasta)     // If aligner is bwa-mem2
        DRAGMAP_HASHTABLE(fasta) // If aligner is dragmap

        bwa         = BWAMEM1_INDEX.out.index.collect()        // path: bwa/*
        bwamem2     = BWAMEM2_INDEX.out.index.collect()        // path: bwamem2/*
        hashtable   = DRAGMAP_HASHTABLE.out.hashmap.collect()    // path: dragmap/*


        versions = versions.mix(BWAMEM1_INDEX.out.versions_bwa)
        versions = versions.mix(BWAMEM2_INDEX.out.versions_bwamem2)
        versions = versions.mix(DRAGMAP_HASHTABLE.out.versions_dragmap)

    } else {

        bwa        = Channel.empty()
        bwamem2    = Channel.empty()
        hashtable  = Channel.empty()
    }

    GATK4_CREATESEQUENCEDICTIONARY(fasta)
    SAMTOOLS_FAIDX(
        fasta.map { meta, fa -> [meta, fa, []] },
        false
    )

    // Canonical nf-core replacement for the deprecated tabix-only wrapper.
    // This normalizes any VCF/VCF.GZ input into a fresh bgzipped VCF plus matching TBI.
    TABIX_DBSNP(
        dbsnp.flatten().map { file -> [[id: file.name.replaceFirst(/\.vcf(\.gz)?$/, '')], file, [], []] },
        'compress',
        true,
        'vcf'
    )
    TABIX_GERMLINE_RESOURCE(
        germline_resource.flatten().map { file -> [[id: file.name.replaceFirst(/\.vcf(\.gz)?$/, '')], file, [], []] },
        'compress',
        true,
        'vcf'
    )
    TABIX_KNOWN_SNPS(
        known_snps.flatten().map { file -> [[id: file.name.replaceFirst(/\.vcf(\.gz)?$/, '')], file, [], []] },
        'compress',
        true,
        'vcf'
    )
    TABIX_KNOWN_INDELS(
        known_indels.flatten().map { file -> [[id: file.name.replaceFirst(/\.vcf(\.gz)?$/, '')], file, [], []] },
        'compress',
        true,
        'vcf'
    )
    TABIX_PON(
        pon.flatten().map { file -> [[id: file.name.replaceFirst(/\.vcf(\.gz)?$/, '')], file, [], []] },
        'compress',
        true,
        'vcf'
    )

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.rna){
        if (params.gtf) {
            if (params.gtf.endsWith('.gz')) {
                GUNZIP_GTF (
                    Channel.fromPath(params.gtf).map{ it -> [[id:it[0].baseName], it] }
                )
                ch_gtf = GUNZIP_GTF.out.gunzip.collect()
                versions = versions.mix(GUNZIP_GTF.out.versions_gunzip)
            } else {
                ch_gtf = Channel.fromPath(params.gtf).collect().map{gtf -> [[id:"gtf"], gtf]}
            }
        } else if (params.gff) {
            if (params.gff.endsWith('.gz')) {
                GUNZIP_GFF (
                    Channel.fromPath(params.gff).map{ it -> [[id:it[0].baseName], it] }
                )
                ch_gff = GUNZIP_GFF.out.gunzip.collect()
                versions = versions.mix(GUNZIP_GFF.out.versions_gunzip)
            } else {
                ch_gff = Channel.fromPath(params.gff).collect().map{gff -> [[id:"gff"], gff]}
            }

            GFFREAD (
                ch_gff,
                fasta
            )
            .gtf
            .set { ch_gtf }

            versions = versions.mix(GFFREAD.out.versions_gffread)
        }

        //
        // Uncompress STAR index or generate from scratch if required
        //
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                UNTAR_STAR_INDEX (
                    Channel.fromPath(params.star_index).map{ it -> [[id:it[0].baseName], it] }
                )
                ch_star_index = UNTAR_STAR_INDEX.out.untar.collect()
                versions   = versions.mix(UNTAR_STAR_INDEX.out.versions_untar)
            } else {
                ch_star_index = Channel.fromPath(params.star_index).collect().map{star_index -> [[id:"star_index"], star_index]}
            }
        }
        else {
            STAR_GENOMEGENERATE ( fasta, ch_gtf )
            ch_star_index = STAR_GENOMEGENERATE.out.index
            versions      = versions.mix(STAR_GENOMEGENERATE.out.versions_star)
        }


        // HISAT2 not necessary if second pass skipped
        if ((params.tools && params.tools.split(',').contains("realignment"))){
            if (params.splicesites) {
                ch_splicesites  = Channel.fromPath(params.splicesites).collect().map{ it -> [ [ id:'null' ], it ]}
            } else{
                HISAT2_EXTRACTSPLICESITES ( ch_gtf )
                ch_splicesites  = HISAT2_EXTRACTSPLICESITES.out.txt
                versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions_hisat2)
            }

            if (params.hisat2_index) {
                ch_hisat2_index  = Channel.fromPath(params.hisat2_index).collect().map{it -> [ [ id:"hisat2_index" ], it ]}
            } else{
                HISAT2_BUILD (
                                fasta,
                                ch_gtf,
                                ch_splicesites
                            )
                ch_hisat2_index = HISAT2_BUILD.out.index
                versions = versions.mix(HISAT2_BUILD.out.versions_hisat2)
            }
        } else {
            ch_hisat2_index = []
            ch_splicesites  = []
        }
    } else {
        ch_star_index   = Channel.empty()
        ch_gtf          = Channel.empty()
        ch_hisat2_index = Channel.empty()
        ch_splicesites  = Channel.value([])
    }


    // Gather versions of all tools used
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions_samtools)
    versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions_gatk4)
    versions = versions.mix(TABIX_DBSNP.out.versions_htslib)
    versions = versions.mix(TABIX_DBSNP.out.versions_xz)
    versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions_htslib)
    versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions_xz)
    versions = versions.mix(TABIX_KNOWN_SNPS.out.versions_htslib)
    versions = versions.mix(TABIX_KNOWN_SNPS.out.versions_xz)
    versions = versions.mix(TABIX_KNOWN_INDELS.out.versions_htslib)
    versions = versions.mix(TABIX_KNOWN_INDELS.out.versions_xz)
    versions = versions.mix(TABIX_PON.out.versions_htslib)
    versions = versions.mix(TABIX_PON.out.versions_xz)

    emit:
        bwa                   = bwa       // path: bwa/*
        bwamem2               = bwamem2       // path: bwamem2/*
        hashtable             = hashtable // path: dragmap/*
        dbsnp                 = TABIX_DBSNP.out.output.map{ meta, gz -> [gz] }.collect()
        dbsnp_tbi             = TABIX_DBSNP.out.index.map{ meta, tbi -> [tbi] }.collect()               // path: dbsnb.vcf.gz.tbi
        dict                  = GATK4_CREATESEQUENCEDICTIONARY.out.dict                               // path: genome.fasta.dict
        fasta                 = fasta
        fasta_fai             = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }                      // path: genome.fasta.fai
        germline_resource     = TABIX_GERMLINE_RESOURCE.out.output.map{ meta, gz -> [gz] }.collect()
        germline_resource_tbi = TABIX_GERMLINE_RESOURCE.out.index.map{ meta, tbi -> [tbi] }.collect()   // path: germline_resource.vcf.gz.tbi
        known_snps            = TABIX_KNOWN_SNPS.out.output.map{ meta, gz -> [gz] }.collect()
        known_snps_tbi        = TABIX_KNOWN_SNPS.out.index.map{ meta, tbi -> [tbi] }.collect()          // path: {known_indels*}.vcf.gz.tbi
        known_indels          = TABIX_KNOWN_INDELS.out.output.map{ meta, gz -> [gz] }.collect()
        known_indels_tbi      = TABIX_KNOWN_INDELS.out.index.map{ meta, tbi -> [tbi] }.collect()        // path: {known_indels*}.vcf.gz.tbi
        pon                   = TABIX_PON.out.output.map{ meta, gz -> [gz] }.collect()
        pon_tbi               = TABIX_PON.out.index.map{ meta, tbi -> [tbi] }.collect()                 // path: pon.vcf.gz.tbi
        star_index            = ch_star_index                                                         // path: star/index/
        gtf                   = ch_gtf                                                                // path: genome.gtf
        hisat2_index          = ch_hisat2_index                                                       // path: hisat2/index/
        splicesites           = ch_splicesites                                                        // path: genome.splicesites.txt
        versions              = versions                                                              // channel: [ versions.yml ]
}
