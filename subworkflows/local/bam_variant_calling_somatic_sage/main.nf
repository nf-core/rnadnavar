//
// SAGE variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { SAGE                                   } from '../../../modules/local/sage/main'
include { GATK4_MERGEVCFS as MERGE_SAGE          } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { HTSLIB_BGZIPTABIX as BGZIPTABIX_VC_SAGE } from '../../../modules/nf-core/htslib/bgziptabix/main'
include { UNZIP as UNZIP_SAGE_ENSEMBL            } from '../../../modules/nf-core/unzip/main'
include { UNTAR as UNTAR_SAGE_ENSEMBL            } from '../../../modules/nf-core/untar/main'

workflow BAM_VARIANT_CALLING_SOMATIC_SAGE {
    take:
    cram      // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi ] manta* are optional
    dict      // channel: [mandatory] [ meta, dict ]
    fasta     // channel: [mandatory] [ meta, fasta ]
    fasta_fai // channel: [mandatory] [ meta, fasta_fai ]
    intervals // channel: [mandatory] [ intervals, num_intervals ] or [ [], 0 ] if no intervals

    main:
    versions = Channel.empty()

    // Unzip/untar resource files if needed
    // prepare vep and sage resource files
    if (!params.sage_ensembl_dir) sage_ensembl = Channel.value([])
    else if (params.sage_ensembl_dir.endsWith(".tar.gz")) {
        UNTAR_SAGE_ENSEMBL(Channel.fromPath(params.sage_ensembl_dir).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        sage_ensembl = UNTAR_SAGE_ENSEMBL.out.untar
        versions = versions.mix(UNTAR_SAGE_ENSEMBL.out.versions_untar)
    } else if (params.sage_ensembl_dir.endsWith(".zip")) {
        UNZIP_SAGE_ENSEMBL(Channel.fromPath(params.sage_ensembl_dir).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
        sage_ensembl = UNZIP_SAGE_ENSEMBL.out.unzipped_archive
        versions = versions.mix(UNZIP_SAGE_ENSEMBL.out.versions_7za)
    } else {
        sage_ensembl = Channel.fromPath(params.sage_ensembl_dir).collect().map{ it -> [ [ id:it[0].baseName ], it ] }
    }

    sage_ensembl.dump(tag:"sage_ensembl")
    // Combine cram and intervals for spread and gather strategy
    cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map and reorganize channel for SAGE module
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, intervls, num_intervals -> [ meta + [ num_intervals:num_intervals ], normal_cram, normal_crai, tumor_cram, tumor_crai, intervls ]}
    SAGE(
        cram_intervals,
        sage_ensembl,
        Channel.fromPath(params.sage_high_confidence).collect().map{ it -> [ [ id:it[0].baseName ], it ] },
        Channel.fromPath(params.sage_actionable_panel).collect().map{ it -> [ [ id:it[0].baseName ], it ] },
        Channel.fromPath(params.sage_known_hotspots).collect().map{ it -> [ [ id:it[0].baseName ], it ] },
        fasta,
        fasta_fai,
        dict)

    // Figuring out if there is one or more vcf(s) from the same sample
    sage_vcf_out = SAGE.out.vcf.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals
    vcf_to_merge = sage_vcf_out.intervals.map{ meta, vcf -> [ groupKey(meta, meta.num_intervals), vcf ]}.groupTuple()
    MERGE_SAGE(vcf_to_merge, dict)
    merged_vcf = MERGE_SAGE.out.vcf.join(MERGE_SAGE.out.tbi, failOnMismatch: true)
    // Only when no_intervals
    BGZIPTABIX_VC_SAGE(
        sage_vcf_out.no_intervals.map{ meta, vcf -> [ meta, vcf, [], [] ] },
        'compress',
        true,
        'vcf'
    )

    // Mix intervals and no_intervals channels together
    vcf = BGZIPTABIX_VC_SAGE.out.output
        .join(BGZIPTABIX_VC_SAGE.out.index, failOnMismatch: true)
        .mix(merged_vcf)
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf, tbi -> [ meta - meta.subMap('normal_id', 'tumor_id','num_intervals') + [ variantcaller:'sage' ], vcf ] }

    versions = versions.mix(MERGE_SAGE.out.versions_gatk4)
    versions = versions.mix(SAGE.out.versions)
    versions = versions.mix(BGZIPTABIX_VC_SAGE.out.versions_htslib)
    versions = versions.mix(BGZIPTABIX_VC_SAGE.out.versions_xz)

    emit:

    vcf
    versions
}
