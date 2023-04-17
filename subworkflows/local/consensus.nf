//
// Make a consensus out of different VCF files
//

include { RUN_CONSENSUS                                    } from '../../modules/local/run_consensus'
include { RUN_CONSENSUS as RUN_CONSENSUS_RESCUE_DNA        } from '../../modules/local/run_consensus'
include { RUN_CONSENSUS as RUN_CONSENSUS_RESCUE_RNA        } from '../../modules/local/run_consensus'
include { PAIR_VARIANT_CALLING_MUTECT2 as GATK_FORCE_CALLS } from './force_mutect_pair_variant_calling'

include { GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING as GATK_TUMOR_ONLY_SOMATIC_VARIANT_CALLING_FORCE } from '../nf-core/gatk4/tumor_only_somatic_variant_calling/main'

workflow CONSENSUS {
    take:
        vcfs
        ch_cram_variant_calling  // from previous variant calling
        dict
        fasta
        fasta_fai
        germline_resource
        germline_resource_tbi
        intervals
        intervals_bed_gz_tbi
        intervals_bed_combined
        pon
        pon_tbi
        previous_vcf_consensus_dna  // results already done to avoid a second run when rna filterig
        previous_vcfs_status_dna  // results already done to avoid a second run when rna filterig
        ch_input_sample
        contamination
        segmentation
        orientation

    main:
        ch_versions          = Channel.empty()
        ch_reports           = Channel.empty()
        // CONSENSUS
        if (params.step in ['mapping', 'markduplicates', 'splitncigar', 'prepare_recalibration', 'recalibrate', 'variant_calling', 'normalize', 'consensus'] ) {
            // set inputs
            if (params.step == 'consensus') {
                vcf_to_consensus = ch_input_sample // TODO: testing
            } else {
                vcf_to_consensus = vcfs
            }
            // Note that DNA can be done from before so we do not need to re-do it
            if (previous_vcf_consensus_dna){ // if the value is given
                // we will only do consensus for RNA again and not DNA again as results are given
                vcf_to_consensus.branch{
                                        dna: it[0].status < 2
                                        rna: it[0].status == 2
                                         }.set{vcf_to_consensus_status}
                vcf_to_consensus = vcf_to_consensus_status.rna  // if second pass we only do CONSENSUS to RNA
                vcf_consensus_status_dna = previous_vcf_consensus_dna   // VCF with consensus calling
                vcfs_status_dna = previous_vcfs_status_dna   // VCF with consensus calling
                vcf_to_consensus.dump(tag:"[STEP5 CONSENSUS] - vcf_to_consensus_2pass --")
                vcfs_status_dna.dump(tag:"[STEP5 CONSENSUS] - vcfs_status_dna_2pass --")
                vcf_consensus_status_dna.dump(tag:"[STEP5 CONSENSUS] - vcf_consensus_status_dna_2pass --")
            } else{ // if the value is null we initialize the channels
                vcf_consensus_dna        = null
                vcfs_status_dna          = null
                vcf_consensus_status_dna = null
            }

            // group vcf per id, sample and status
            vcf_to_consensus.dump(tag:"[STEP5 CONSENSUS] - vcf_to_consensus0")
//            ncallers = Channel.from(list.count('sage') + list.count('strelka') + list.count('mutect2') + list.count("consensus"))

            vcf_to_consensus = vcf_to_consensus.map{ meta, vcf ->
                                        def toolsllist = params.tools.split(',')
                                        def ncallers = toolsllist.count('sage') +
                                                       toolsllist.count('strelka') +
                                                       toolsllist.count('mutect2') +
                                                       toolsllist.count("consensus") +
                                                       toolsllist.count("freebayes")
                                        [groupKey([id:meta.id,
                                             patient:meta.patient,
                                             status:meta.status,
                                             ncallers:ncallers
                                             ], ncallers), vcf, meta.variantcaller]}
                                     .groupTuple() // makes the whole pipeline wait for all processes to finish

            vcf_to_consensus.dump(tag:"[STEP5 CONSENSUS] - vcf_to_consensus1")
            RUN_CONSENSUS ( vcf_to_consensus )
            ch_versions   = RUN_CONSENSUS.out.versions
            consensus_vcf = RUN_CONSENSUS.out.vcf  // 1 consensus_vcf from all callers

            // CONSENSUS 2nd pass if DNA and RNA - RESCUE STEP TODO: not sure if it works when only 1
            // 1. Separate DNA from RNA
            //  For consensus vcf
            consensus_vcf.branch{
                dna: it[0].status < 2
                rna: it[0].status == 2
            }.set{vcf_consensus_status}
            // For input vcfs -- consensus vcf needs to be added
            vcf_to_consensus.branch{
                dna: it[0].status < 2
                rna: it[0].status == 2
            }.set{vcfs_status}

            vcfs_status.dna.dump(tag:"[STEP5 CONSENSUS] - vcfs_status.dna")
            vcf_consensus_status.dna.dump(tag:"[STEP5 CONSENSUS] - vcf_consensus_status.dna")

            if(!previous_vcf_consensus_dna){
                vcf_consensus_dna        = vcf_consensus_status.dna
                vcfs_status_dna          = vcfs_status.dna
                vcf_consensus_status_dna = vcf_consensus_status.dna
            }

            vcf_consensus_status_dna.dump(tag:"[STEP5 CONSENSUS] vcf_consensus_status_dna")
            vcfs_status_dna.dump(tag:"[STEP5 CONSENSUS] vcfs_status_dna")

        if (params.skip_tools && params.skip_tools.split(',').contains('rescue')) {
            // no rescue step
            // TODO TEST this code
            vcf_consensus_status_dna = vcf_consensus_status.dna
            vcf_consensus_status_rna = vcf_consensus_status.rna
        } else {
            // RESCUE STEP
            // cross dna / rna
            // VCF from consensus
            vcf_consensus_status_dna_to_cross = vcf_consensus_status_dna.map{
                                            meta, vcf, caller ->
                                            [meta.patient, meta, [vcf], caller] }

            vcf_consensus_status_rna_to_cross = vcf_consensus_status.rna.map{
                                            meta, vcf, caller ->
                                            [meta.patient, meta, [vcf], caller] }
            // VCFs from variant calling
            vcfs_status_dna_to_cross = vcfs_status_dna.map{
                                            meta, vcfs, callers ->
                                            [meta.patient, meta, vcfs, callers] }

            vcfs_status_rna_to_cross = vcfs_status.rna.map{
                                            meta, vcfs, callers ->
                                            [meta.patient, meta, vcfs, callers] }

            vcf_consensus_status_dna_to_cross.dump(tag:"[STEP5 CONSENSUS] vcf_consensus_status_dna_to_cross")
            vcf_consensus_status_rna_to_cross.dump(tag:"[STEP5 CONSENSUS] vcf_consensus_status_rna_to_cross")

            // cross results keeping metadata
            vcfs_dna_crossed_with_rna_rescue = vcfs_status_dna_to_cross
                                               .cross(vcf_consensus_status_rna_to_cross)
                                               .map { dna, rna ->
                                                def meta = [:]
                                                meta.patient = dna[0]
                                                meta.dna_id  = dna[1].id
                                                meta.rna_id  = rna[1].id
                                                meta.status  = dna[1].status
                                                meta.id      = "${meta.dna_id}_with_${meta.rna_id}".toString()
                                                [meta, dna[2] + rna[2], dna[3] + rna[3]]
                                            }
            vcfs_rna_crossed_with_dna_rescue = vcfs_status_rna_to_cross
                                               .cross(vcf_consensus_status_dna_to_cross)
                                               .map { rna, dna ->
                                                def meta = [:]
                                                meta.patient = rna[0]
                                                meta.rna_id  = rna[1].id
                                                meta.dna_id  = dna[1].id
                                                meta.status  = rna[1].status
                                                meta.id      = "${meta.rna_id}_with_${meta.dna_id}".toString()
                                                [meta, rna[2] + dna[2], rna[3] + dna[3]]
                                            }


            // TODO: unify?
            vcfs_dna_crossed_with_rna_rescue.dump(tag:"[STEP5 CONSENSUS] vcfs_dna_crossed_with_rna_rescue")
            vcfs_rna_crossed_with_dna_rescue.dump(tag:"[STEP5 CONSENSUS] vcfs_rna_crossed_with_dna_rescue")
            RUN_CONSENSUS_RESCUE_DNA ( vcfs_dna_crossed_with_rna_rescue )
            RUN_CONSENSUS_RESCUE_RNA ( vcfs_rna_crossed_with_dna_rescue )

            vcf_consensus_dna = RUN_CONSENSUS_RESCUE_DNA.out.vcf
            vcf_consensus_rna = RUN_CONSENSUS_RESCUE_RNA.out.vcf

        }

        }

    emit:
        vcf_consensus_dna   = vcf_consensus_dna
        vcfs_status_dna     = vcfs_status_dna
        maf                 = vcf_consensus_dna.mix(vcf_consensus_rna) // channel: [ [meta], vcf ]
        reports             = ch_reports // channel: [ versions.yml ]
        versions            = ch_versions // channel: [ versions.yml ]
}