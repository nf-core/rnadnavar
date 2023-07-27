//
// ANNOTATION
//

include { ENSEMBLVEP       } from '../../modules/nf-core/modules/ensemblvep/main'
include { VCF2MAF          } from '../../modules/local/vcf2maf/vcf2maf/main'


workflow ANNOTATE {
    take:
        tools
        vcf_to_annotate          // channel: [ val(meta), vcf ]
        fasta
        ch_input_sample

    main:
        ch_reports   = Channel.empty()
        ch_vcf_ann   = Channel.empty()
        ch_maf_ann   = Channel.empty()
        ch_tab_ann   = Channel.empty()
        ch_json_ann  = Channel.empty()
        ch_versions  = Channel.empty()

        // Initialize value channels based on params, defined in the params.genomes[params.genome] scope
        vep_cache_version  = params.vep_cache_version  ?: Channel.empty()
        vep_genome         = params.vep_genome         ?: Channel.empty()
        vep_species        = params.vep_species        ?: Channel.empty()
        vep_cache          = params.vep_cache          ? Channel.fromPath(params.vep_cache).collect()                : []

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

        vep_fasta = (params.vep_include_fasta) ? fasta : []

        if (params.step == 'annotate') {
            vcf_to_annotate = ch_input_sample
        }
        if (tools.split(',').contains('vep')) {
            ENSEMBLVEP(vcf_to_annotate, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)
            ch_versions = ch_versions.mix(ENSEMBLVEP.out.versions.first())
            ch_reports   = ch_reports.mix(ENSEMBLVEP.out.report)

            ch_vcf_ann   = ch_vcf_ann.mix(ENSEMBLVEP.out.vcf)
            ch_tab_ann   = ch_vcf_ann.mix(ENSEMBLVEP.out.tab)
            ch_json_ann  = ch_vcf_ann.mix(ENSEMBLVEP.out.json)
        }
        if (tools.split(',').contains('vcf2maf')) {
            VCF2MAF(ch_vcf_ann,
                    fasta)
            ch_maf_ann = VCF2MAF.out.maf
        }


    emit:
        maf_ann   = ch_maf_ann      // channel: [ val(meta), maf ]
        vcf_ann   = ch_vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
        tab_ann   = ch_tab_ann
        json_ann  = ch_json_ann
        reports   = ch_reports      //    path: *.html
        versions  = ch_versions     //    path: versions.yml
}
