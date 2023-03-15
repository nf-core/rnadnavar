//
// Make a consensus out of different VCF files
//

include { RUN_CONSENSUS                                    } from '../../modules/local/run_consensus'
include { RUN_CONSENSUS as RUN_CONSENSUS_RESCUE_DNA        } from '../../modules/local/run_consensus'
include { RUN_CONSENSUS as RUN_CONSENSUS_RESCUE_RNA        } from '../../modules/local/run_consensus'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_DNA         } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_RNA         } from '../../modules/nf-core/modules/tabix/bgziptabix/main'
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

            vcf_to_consensus = vcf_to_consensus.map{ meta, vcf->
                                        [[id:meta.id,
                                         patient:meta.patient,
                                         status:meta.status,
                                         ncallers:
                                         ], vcf, meta.variantcaller]
                                        }.groupTuple(size:3) // makes the whole pipeline wait for all processes to finish

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

            contamination.branch{ dna: it[0].status < 2
                                   rna: it[0].status == 2}
                          .set{contamination_status}
            segmentation.branch{ dna: it[0].status < 2
                                   rna: it[0].status == 2}
                          .set{segmentation_status}
            orientation.branch{ dna: it[0].status < 2
                                   rna: it[0].status == 2}
                          .set{orientation_status}

//            contamination = contamination.map{meta, table ->[meta.patient, meta, table] }
//            segmentation  = segmentation.map{meta, table ->[meta.patient, meta, table] }
//            orientation   = orientation.map{meta, table ->[meta.patient, meta, table] }
            contamination_status_dna = contamination_status.dna.map{meta, table -> [meta.patient, meta, table]}
            contamination_status_rna = contamination_status.rna.map{meta, table -> [meta.patient, meta, table]}
            segmentation_status_dna  = segmentation_status.dna.map{meta, table -> [meta.patient, meta, table]}
            segmentation_status_rna  = segmentation_status.rna.map{meta, table -> [meta.patient, meta, table]}
            orientation_status_dna   = orientation_status.dna.map{meta, table -> [meta.patient, meta, table]}
            orientation_status_rna   = orientation_status.rna.map{meta, table -> [meta.patient, meta, table]}
            contamination_status_dna.dump(tag:"contamination_status_dna")
            // TODO check the [] from before!
//            rna_contamination_with_dna_rescue = contamination_status_rna
//                                               .cross(vcf_consensus_status_dna_to_cross)
//                                               .map { rna, dna ->
//                                                def meta = [:]
//                                                // add alleles
//                                                meta.patient = rna[0]
//                                                meta.rna_id  = rna[1].id
//                                                meta.dna_id  = dna[1].id
//                                                meta.status  = rna[1].status
//                                                meta.id      = "${meta.rna_id}_with_${meta.dna_id}".toString()
//                                                [ rna[0] + [], rna[2]]
//                                            }
//            dna_contamination_with_rna_rescue = contamination_status_dna
//                                               .cross(vcf_consensus_status_rna_to_cross)
//                                               .map { dna, rna ->
//                                               def meta = [:]
//                                                meta.patient = dna[0]
//                                                meta.dna_id  = dna[1].id
//                                                meta.rna_id  = rna[1].id
//                                                meta.status  = dna[1].status
//                                                meta.id      = "${meta.dna_id}_with_${meta.rna_id}".toString()
//                                                [meta, dna[2]]
//                                            }
//
//            rna_segmentation_with_dna_rescue = segmentation_status_rna
//                                               .cross(vcf_consensus_status_dna_to_cross)
//                                               .map { rna, dna ->
//                                                def meta = [:]
//                                                meta.patient = rna[0]
//                                                meta.rna_id  = rna[1].id
//                                                meta.dna_id  = dna[1].id
//                                                meta.status  = rna[1].status
//                                                meta.id      = "${meta.rna_id}_with_${meta.dna_id}".toString()
//                                                [meta, rna[2]]
//                                            }
//            dna_segmentation_with_rna_rescue = segmentation_status_dna
//                                                .cross(vcf_consensus_status_rna_to_cross)
//                                                .map { dna, rna ->
//                                                def meta = [:]
//                                                meta.patient = dna[0]
//                                                meta.dna_id  = dna[1].id
//                                                meta.rna_id  = rna[1].id
//                                                meta.status  = dna[1].status
//                                                meta.id      = "${meta.dna_id}_with_${meta.rna_id}".toString()
//                                                [meta, dna[2]]
//                                            }
//            segmentation_with_rescue = dna_segmentation_with_rna_rescue.mix(rna_segmentation_with_dna_rescue)
//
//            rna_orientation_with_dna_rescue = orientation_status_rna
//                                               .cross(vcf_consensus_status_dna_to_cross)
//                                               .map { rna, dna ->
//                                                def meta = [:]
//                                                meta.patient = rna[0]
//                                                meta.rna_id  = rna[1].id
//                                                meta.dna_id  = dna[1].id
//                                                meta.status  = rna[1].status
//                                                meta.id      = "${meta.rna_id}_with_${meta.dna_id}".toString()
//                                                [meta, rna[2]]
//                                            }
//            dna_orientation_with_rna_rescue = orientation_status_dna
//                                               .cross(vcf_consensus_status_rna_to_cross)
//                                               .map { dna, rna ->
//                                               def meta = [:]
//                                                meta.patient = dna[0]
//                                                meta.dna_id  = dna[1].id
//                                                meta.rna_id  = rna[1].id
//                                                meta.status  = dna[1].status
//                                                meta.id      = "${meta.dna_id}_with_${meta.rna_id}".toString()
//                                                [meta, dna[2]]
//                                            }
//            orientation_with_rescue = dna_orientation_with_rna_rescue.mix(rna_orientation_with_dna_rescue)
//
//            contamination_status.rna.dump(tag:"contamination_status")
//            contamination_status.rna.dump(tag:"contamination_status")
//            rna_contamination_with_dna_rescue.dump(tag:"TESTING")

            // TODO: unify?
//            vcfs_crossed_for_rescue = vcfs_dna_crossed_with_rna_rescue.mix(vcfs_rna_crossed_with_dna_rescue)
            vcfs_dna_crossed_with_rna_rescue.dump(tag:"[STEP5 CONSENSUS] vcfs_dna_crossed_with_rna_rescue")
            vcfs_rna_crossed_with_dna_rescue.dump(tag:"[STEP5 CONSENSUS] vcfs_rna_crossed_with_dna_rescue")
            RUN_CONSENSUS_RESCUE_DNA ( vcfs_dna_crossed_with_rna_rescue )
            RUN_CONSENSUS_RESCUE_RNA ( vcfs_rna_crossed_with_dna_rescue )

            vcf_consensus_dna = RUN_CONSENSUS_RESCUE_DNA.out.vcf
            vcf_consensus_rna = RUN_CONSENSUS_RESCUE_RNA.out.vcf

//            TABIX_BGZIPTABIX_DNA (vcf_consensus_dna)
//            TABIX_BGZIPTABIX_RNA (vcf_consensus_rna)
//            vcf_consensus_dna_tbi = TABIX_BGZIPTABIX_DNA.out.gz_tbi
//            vcf_consensus_rna_tbi = TABIX_BGZIPTABIX_RNA.out.gz_tbi
//
//            // END OF RESCUE STEP? Can we use just the consensus without the rescue?
//
//            // Now we have the vfc consensus, however it does not have any annotation.
//            // We force-run GATK variant calling again on these variants only to unify format
//            ch_cram_variant_calling.branch{
//                                           dna: it[0].status < 2
//                                           rna: it[0].status == 2
//                                           }
//                                    .set{ch_cram_variant_calling_status}
//            // create alleles channel with vcf from consensus per patient
//            alleles_dna = vcf_consensus_dna_tbi.map{meta, vcf, tbi -> [meta.patient, vcf]}
//            alleles_rna = vcf_consensus_rna_tbi.map{meta, vcf, tbi -> [meta.patient, vcf]}
//            alleles_dna.dump(tag:"[STEP5 CONSENSUS] alleles_dna")
//            alleles_rna.dump(tag:"[STEP5 CONSENSUS] alleles_rna")
//
//            ch_cram_variant_calling_dna_pair_con = ch_cram_variant_calling_status.dna
//                .map{meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
//                            [meta.patient, meta, normal_cram, normal_crai, tumor_cram, tumor_crai]}
//                .cross(alleles_rna)
//                .map{ crams, alleles ->
//                           def meta = [:]
//                            meta.patient    = crams[0]
//                            meta.normal_id  = crams[1].normal_id
//                            meta.tumor_id   = crams[1].tumor_id
//                            meta.status     = crams[1].status
//                            meta.id         = crams[1].id
//                            meta.alleles    = alleles[1]
//                            [meta, crams[2], crams[3], crams[4], crams[5]]
//                    }
//            ch_cram_variant_calling_rna_pair_con = ch_cram_variant_calling_status.rna
//                .map{meta, normal_cram, normal_crai, tumor_cram, tumor_crai ->
//                            [meta.patient, meta, normal_cram, normal_crai, tumor_cram, tumor_crai]}
//                .cross(alleles_dna)
//                .map{ crams, alleles ->
//                           def meta = [:]
//                            meta.patient    = crams[0]
//                            meta.normal_id  = crams[1].normal_id
//                            meta.tumor_id   = crams[1].tumor_id
//                            meta.status     = crams[1].status
//                            meta.id         = crams[1].id
//                            meta.alleles    = alleles[1]
//                            [meta, crams[2], crams[3], crams[4], crams[5]]
//                    }
//            rna_contamination_with_dna_rescue = contamination_status.rna
//                .map{meta, table ->
//                            [meta.patient, meta, table]}
//                .cross(alleles_dna)
//                .map{ table, alleles ->
//                           def meta = [:]
//                            meta.patient    = table[0]
//                            meta.normal_id  = table[1].normal_id
//                            meta.tumor_id   = table[1].tumor_id
//                            meta.status     = table[1].status
//                            meta.id         = table[1].id
//                            meta.alleles    = alleles[1]
//                            [meta, table[2]]
//                    }
//            dna_contamination_with_rna_rescue = contamination_status.dna
//                .map{meta, table ->
//                            [meta.patient, meta, table]}
//                .cross(alleles_rna)
//                .map{ table, alleles ->
//                           def meta = [:]
//                            meta.patient    = table[0]
//                            meta.normal_id  = table[1].normal_id
//                            meta.tumor_id   = table[1].tumor_id
//                            meta.status     = table[1].status
//                            meta.id         = table[1].id
//                            meta.alleles    = alleles[1]
//                            [meta, table[2]]
//                    }
//            contamination_with_rescue = dna_contamination_with_rna_rescue.mix(rna_contamination_with_dna_rescue)
//
//            // ABOVE is the solution to che contamination problems we need to combine crams info with alleles now plus the contamination
//            rna_segmentation_with_dna_rescue = segmentation_status.rna
//                .map{meta, table ->
//                            [meta.patient, meta, table]}
//                .cross(alleles_dna)
//                .map{ table, alleles ->
//                           def meta = [:]
//                            meta.patient    = table[0]
//                            meta.normal_id  = table[1].normal_id
//                            meta.tumor_id   = table[1].tumor_id
//                            meta.status     = table[1].status
//                            meta.id         = table[1].id
//                            meta.alleles    = alleles[1]
//                            [meta, table[2]]
//                    }
//            dna_segmentation_with_rna_rescue = segmentation_status.dna
//                .map{meta, table ->
//                            [meta.patient, meta, table]}
//                .cross(alleles_rna)
//                .map{ table, alleles ->
//                           def meta = [:]
//                            meta.patient    = table[0]
//                            meta.normal_id  = table[1].normal_id
//                            meta.tumor_id   = table[1].tumor_id
//                            meta.status     = table[1].status
//                            meta.id         = table[1].id
//                            meta.alleles    = alleles[1]
//                            [meta, table[2]]
//                    }
//            segmentation_with_rescue = dna_segmentation_with_rna_rescue.mix(rna_contamination_with_dna_rescue)
//
//            rna_orientation_with_dna_rescue = orientation_status.rna
//                .map{meta, table ->
//                            [meta.patient, meta, table]}
//                .cross(alleles_dna)
//                .map{ table, alleles ->
//                           def meta = [:]
//                            meta.patient    = table[0]
//                            meta.normal_id  = table[1].normal_id
//                            meta.tumor_id   = table[1].tumor_id
//                            meta.status     = table[1].status
//                            meta.id         = table[1].id
//                            meta.alleles    = alleles[1]
//                            [meta, table[2]]
//                    }
//            dna_orientation_with_rna_rescue = orientation_status.dna
//                .map{meta, table ->
//                            [meta.patient, meta, table]}
//                .cross(alleles_rna)
//                .map{ table, alleles ->
//                           def meta = [:]
//                            meta.id         = table[1].id
//                            meta.normal_id  = table[1].normal_id
//                            meta.patient    = table[0]
//                            meta.num_intervals = table[1].num_intervals
//                            meta.tumor_id   = table[1].tumor_id
//                            meta.status     = table[1].status
//                            meta.alleles    = alleles[1]
//                            [meta, table[2]]
//                    }
//            orientation_with_rescue = dna_orientation_with_rna_rescue.mix(rna_contamination_with_dna_rescue)
//
//            contamination_with_rescue.dump(tag:"[STEP5 CONSENSUS] contamination_with_rescue")
//            // cram from used in variant calling with alleles of the consensus will be used to force calls again
//            ch_cram_variant_calling_pairs_con = ch_cram_variant_calling_dna_pair_con.mix(ch_cram_variant_calling_rna_pair_con)
//            // force calls to get consensus with same annotation
//            GATK_FORCE_CALLS(
//                ch_cram_variant_calling_pairs_con,
//                dict,
//                fasta,
//                fasta_fai,
//                germline_resource,
//                germline_resource_tbi,
//                intervals,
//                intervals_bed_gz_tbi,
//                intervals_bed_combined,
//                pon,
//                pon_tbi,
//                params.skip_tools,
//                contamination_with_rescue,
//                segmentation_with_rescue,
//                orientation_with_rescue
//            )
//            ch_versions = ch_versions.mix(GATK_FORCE_CALLS.out.versions)
//            ch_reports  = ch_reports.mix(GATK_FORCE_CALLS.out.reports)
//            GATK_FORCE_CALLS.out.vcf.dump(tag:"[STEP5 CONSENSUS] forced_vcf")
        }

    emit:
        vcf_consensus_dna   = vcf_consensus_dna
        vcfs_status_dna     = vcfs_status_dna
        maf                 = vcf_consensus_dna.mix(vcf_consensus_rna) // channel: [ [meta], vcf ]
        reports             = ch_reports // channel: [ versions.yml ]
        versions            = ch_versions // channel: [ versions.yml ]
}