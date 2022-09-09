//
// Make a consensus out of different VCF files
//

include { RUN_CONSENSUS } from '../../modules/local/run_consensus'
include { RUN_CONSENSUS as RUN_CONSENSUS_RESCUE_DNA } from '../../modules/local/run_consensus'
include { RUN_CONSENSUS as RUN_CONSENSUS_RESCUE_RNA } from '../../modules/local/run_consensus'


workflow CONSENSUS {
    take:
        vcfs

    main:
        // CONSENSUS
        RUN_CONSENSUS ( vcfs )
        ch_versions = RUN_CONSENSUS.out.versions

        // CONSENSUS 2nd pass if DNA and RNA - RESCUE STEP TODO: not sure if it works when only 1
        // 1. Separate DNA from RNA
        //  For consensus vcf
        RUN_CONSENSUS.out.vcf.branch{
            dna: it[0].status < 2
            rna: it[0].status == 2
        }.set{vcf_consensus_status}
        // For input vcfs -- consensus vcf needs to be added
        vcfs.branch{
            dna: it[0].status < 2
            rna: it[0].status == 2
        }.set{vcf_status}


        // cross dna / rna
        // VCF from consensus
        vcf_consensus_status_dna_to_cross = vcf_consensus_status.dna.map{
                                        meta, vcf, caller ->
                                        [meta.patient, meta, [vcf], caller] }

        vcf_consensus_status_rna_to_cross = vcf_consensus_status.rna.map{
                                        meta, vcf, caller ->
                                        [meta.patient, meta, [vcf], caller] }
        // VCFs from variant calling
        vcfs_status_dna_to_cross = vcf_status.dna.map{
                                        meta, vcfs, callers ->
                                        [meta.patient, meta, vcfs, callers] }

        vcfs_status_rna_to_cross = vcf_status.rna.map{
                                        meta, vcfs, callers ->
                                        [meta.patient, meta, vcfs, callers] }
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

        RUN_CONSENSUS_RESCUE_DNA ( vcfs_dna_crossed_with_rna_rescue )
        RUN_CONSENSUS_RESCUE_RNA ( vcfs_rna_crossed_with_dna_rescue )


   //     consensus = RUN_CONSENSUS_RESCUE_DNA.out.vcf.mix()
    //    rna_consensus = RUN_CONSENSUS_RESCUE_RNA.out.vcf

    //    dna_consensus.dump(tag:'dna_consensus')
    //    dna_consensus.join(original_meta).dump(tag:'dna_consensus')
    // RESCUE
//    CONSENSUS with RNA and vice versa then remove variants and annotate



    emit:
        vcf    = RUN_CONSENSUS.out.vcf // channel: [ [meta], vcf ]
        vcf_separate    = RUN_CONSENSUS.out.vcf_separate // channel: [ [meta], vcf ]
        versions    = ch_versions // channel: [ versions.yml ]
}