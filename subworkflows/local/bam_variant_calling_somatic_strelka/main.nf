//
// STRELKA2 tumor-normal variant calling
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run

include { GATK4_MERGEVCFS as MERGE_STRELKA_INDELS } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA_SNVS   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGE_STRELKA        } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { STRELKA_SOMATIC                         } from '../../../modules/nf-core/strelka/somatic/main'

workflow BAM_VARIANT_CALLING_SOMATIC_STRELKA {
    take:
    cram          // channel: [mandatory] [ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi ] manta* are optional
    dict          // channel: [mandatory] [ meta, dict ]
    fasta         // channel: [mandatory] [ meta, fasta ]
    fasta_fai     // channel: [mandatory] fasta FAI path
    intervals     // channel: [mandatory] [ interval.bed.gz, interval.bed.gz.tbi, num_intervals ] or [ [], [], 0 ] if no intervals
    no_intervals  // true/false

    main:
    versions = Channel.empty()

    // Combine cram and intervals for spread and gather strategy
    if (no_intervals){
        cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, intervals_gz_tbi, num_intervals -> [ meta + [ num_intervals:0 ], normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, [], [] ] }
    } else{
        cram_intervals = cram.combine(intervals)
        // Move num_intervals to meta map
        .map{ meta, normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, intervals_gz_tbi, num_intervals -> [ meta + [ num_intervals:num_intervals ], normal_cram, normal_crai, tumor_cram, tumor_crai, manta_vcf, manta_tbi, intervals_gz_tbi[0], intervals_gz_tbi[1] ] }
    }
    STRELKA_SOMATIC(cram_intervals, fasta.map { _meta, fa -> fa }, fasta_fai)

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_indels = STRELKA_SOMATIC.out.vcf_indels.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Figuring out if there is one or more vcf(s) from the same sample
    vcf_snvs = STRELKA_SOMATIC.out.vcf_snvs.branch{
        // Use meta.num_intervals to asses number of intervals
        intervals:    it[0].num_intervals > 1
        no_intervals: it[0].num_intervals <= 1
    }

    // Only when using intervals TODO: check that commenting map does not break things
    vcfs_to_merge = STRELKA_SOMATIC.out.vcf_indels
                    .mix(STRELKA_SOMATIC.out.vcf_snvs)
                    .map{ meta, vcf -> [ groupKey(meta, meta.num_intervals * 2), vcf ]}
                    .groupTuple()

    // Mix intervals and no_intervals channels together
    // Merge SNVs and indels
    MERGE_STRELKA(vcfs_to_merge, dict)

    // Mix intervals and no_intervals channels together
    vcf = MERGE_STRELKA.out.vcf
        // add variantcaller to meta map and remove no longer necessary field: num_intervals
        .map{ meta, vcf ->
                [ meta.subMap('id', 'patient', 'status', 'variantcaller', 'file_name', 'data_type') + [ variantcaller:'strelka' ], vcf ] }


//    versions = versions.mix(MERGE_STRELKA_SNVS.out.versions)
//    versions = versions.mix(MERGE_STRELKA_INDELS.out.versions)
    versions = versions.mix(MERGE_STRELKA.out.versions_gatk4)
    versions = versions.mix(STRELKA_SOMATIC.out.versions)

    emit:
    vcf      = vcf

    versions
}
