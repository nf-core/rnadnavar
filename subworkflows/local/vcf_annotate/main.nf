//
// ANNOTATION
//

include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep/main'
include { CHANNEL_VARIANT_CALLING_CREATE_CSV as CHANNEL_ANNOTATE_CREATE_CSV  } from '../channel_variant_calling_create_csv/main'

workflow VCF_ANNOTATE {
    take:
    vcf          // channel: [ val(meta), vcf ]
    fasta
    input_sample
    second_run

    main:
    reports  = Channel.empty()
    vcf_ann  = Channel.empty()
    tab_ann  = Channel.empty()
    json_ann = Channel.empty()
    versions = Channel.empty()

    if (params.step == 'annotate') vcf = input_sample

    if (params.tools && params.tools.split(',').contains('vep') || second_run) {

        if (params.tools && params.tools.split(',').contains('vep') || second_run) {
            fasta = (params.vep_include_fasta) ? fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] } : [[id: 'null'], []]
            vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
            vep_genome         = params.vep_genome         ?: Channel.empty()
            vep_species        = params.vep_species        ?: Channel.empty()
            vep_cache          = params.vep_cache    ? params.use_annotation_cache_keys ? Channel.fromPath("${params.vep_cache}/${params.vep_cache_version}_${params.vep_genome}").collect() : Channel.fromPath(params.vep_cache).collect()    : []

            vep_extra_files = []
            if (params.dbnsfp && params.dbnsfp_tbi) {
                vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
                vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
            }

            if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
                vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
                vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
                vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
                vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
            }

            vcf_for_vep = vcf.map{ meta, vcf -> [ meta, vcf, [] ] }
            VCF_ANNOTATE_ENSEMBLVEP(vcf_for_vep, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

            reports  = reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
            vcf_ann  = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf_tbi).map{meta, vcf, tbi -> [meta +[data_type:"vcf"], vcf, tbi]}
            tab_ann  = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
            json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)
            versions = versions.mix(VCF_ANNOTATE_ENSEMBLVEP.out.versions)
            CHANNEL_ANNOTATE_CREATE_CSV(vcf_ann.map{meta, vcf, tbi -> [meta, vcf]}, "annotated")

        }
    } else{

        vcf_ann = vcf.map{metaVCF -> [metaVCF[0] + [data_type:"vcf"], metaVCF[1]]}

    }

    emit:
    vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
    reports      //    path: *.html
    versions     //    path: versions.yml
}
