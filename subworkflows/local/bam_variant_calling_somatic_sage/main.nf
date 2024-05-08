//
// SAGE variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { BCFTOOLS_SORT                         } from '../../../modules/nf-core/bcftools/sort/main'
include { SAGE                                  } from '../../../modules/local/sage/main'
include { GATK4_MERGEVCFS as MERGE_SAGE         } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_TABIX     as TABIX_VC_SAGE      } from '../../../modules/nf-core/tabix/tabix/main'
include { UNZIP as UNZIP_SAGE_ENSEMBL           } from '../../../modules/nf-core/unzip/main'
include { UNTAR as UNTAR_SAGE_ENSEMBL           } from '../../../modules/nf-core/untar/main'

workflow BAM_VARIANT_CALLING_SOMATIC_SAGE {
    take:
    cram      // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi ] manta* are optional
    dict      // channel: [mandatory] [ meta, dict ]
    fasta     // channel: [mandatory] [ fasta ]
    fasta_fai // channel: [mandatory] [ fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Unzip/untar resource files if needed
    // prepare vep and sage resource files
    if (!params.sage_ensembl_dir) sage_ensembl = Channel.value([])
    else if (params.sage_ensembl_dir.endsWith(".tar.gz")) {
        UNTAR_SAGE_ENSEMBL(Channel.fromPath(params.sage_ensembl_dir).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        sage_ensembl = UNTAR_SAGE_ENSEMBL.out.untar
        versions = versions.mix(UNTAR_SAGE_ENSEMBL.out.versions)
    } else if (params.sage_ensembl_dir.endsWith(".zip")) {
        UNZIP_SAGE_ENSEMBL(Channel.fromPath(params.sage_ensembl_dir).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        sage_ensembl = UNZIP_SAGE_ENSEMBL.out.unzipped_archive
        versions = versions.mix(UNZIP_SAGE_ENSEMBL.out.versions)
    } else {
        sage_ensembl = Channel.fromPath(params.sage_ensembl_dir).collect().map{ it -> [ [ id:it[0].baseName ], it ] }
    }

    sage_ensembl.dump(tag:"sage_ensembl")
    // sage_resources.dump(tag:"sage_resources")
    // Combine cram and intervals for spread and gather strategy
    // sage_resources = Channel.value([])
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for SAGE module
        .map{ meta, cram1, crai1, cram2, crai2, intervals, num_intervals -> [ meta + [ num_intervals:num_intervals ], cram1, crai1, cram2, crai2, intervals ]}
    SAGE(
        cram_intervals,
        sage_ensembl,
        Channel.fromPath(params.sage_highconfidence).collect().map{ it -> [ [ id:it[0].baseName ], it ] },
        Channel.fromPath(params.sage_actionablepanel).collect().map{ it -> [ [ id:it[0].baseName ], it ] },
        Channel.fromPath(params.sage_knownhotspots).collect().map{ it -> [ [ id:it[0].baseName ], it ] },
        fasta,
        fasta_fai,
        dict)

    BCFTOOLS_SORT(SAGE.out.vcf)

    // Figuring out if there is one or more vcf(s) from the same sample
    bcftools_vcf_out = BCFTOOLS_SORT.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = bcftools_vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    MERGE_SAGE(vcf_to_merge, dict)

    // Only when no_intervals
    TABIX_VC_SAGE(bcftools_vcf_out.no_intervals)

    // Mix intervals and no_intervals channels together
    vcf = MERGE_SAGE.out.vcf.mix(bcftools_vcf_out.no_intervals)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf -> [ meta - meta.subMap('normal_id', 'tumor_id','num_intervals') + [ variantcaller:'sage' ], vcf ] }

    versions = versions.mix(BCFTOOLS_SORT.out.versions)
    versions = versions.mix(MERGE_SAGE.out.versions)
    versions = versions.mix(SAGE.out.versions)
    versions = versions.mix(TABIX_VC_SAGE.out.versions)

    emit:
    vcf

    versions
}
