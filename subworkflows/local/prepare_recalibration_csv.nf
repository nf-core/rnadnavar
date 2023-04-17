//
// PREPARE_RECALIBRATION_CSV
//

workflow PREPARE_RECALIBRATION_CSV {
    take:
        cram_table_bqsr // channel: [mandatory] meta, cram, crai, table
        skip_tools

    main:
        // Creating csv files to restart from this step
        if (!(skip_tools && (skip_tools.split(',').contains('markduplicates')))) {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, cram, crai, table ->
                patient = meta.patient
                id = meta.id
                lane    = meta.lane
                sample  = meta.sample
                status  = meta.status
                suffix_aligned = params.save_output_as_bam ? "bam" : "cram"
                suffix_index   = params.save_output_as_bam ? "bam.bai" : "cram.crai"
                cram = status <2 ? "${params.outdir}/preprocessing/markduplicates/${patient}/${id}/${cram.baseName}.${suffix_aligned}" : "${params.outdir}/preprocessing/splitncigar/${patient}/${id}/${cram.baseName}.${suffix_aligned}"
                crai = status <2 ? "${params.outdir}/preprocessing/markduplicates/${patient}/${id}/${crai.baseName.minus(".cram")}.${suffix_index}" : "${params.outdir}/preprocessing/splitncigar/${patient}/${id}/${crai.baseName.minus(".cram")}.${suffix_index}"
                table = "${params.outdir}/preprocessing/recal_table/${patient}/${id}/${id}.recal.table"
                ["markduplicates.csv", "patient,status,sample,lane,fastq_1,fastq_2,bam,bai,cram,crai,table,vcf\n${patient},${status},${sample},${lane},,,,,${cram},${crai},${table},\n"]
            }
        } else {
            cram_table_bqsr.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, cram, crai, table ->
                patient = meta.patient
                id = meta.id
                lane    = meta.lane
                sample  = meta.sample
                status  = meta.status
                suffix_aligned = params.save_output_as_bam ? "bam" : "cram"
                suffix_index   = params.save_output_as_bam ? "bam.bai" : "cram.crai"
                cram = status <2 ? "${params.outdir}/preprocessing/mapped/${patient}/${id}/${cram.baseName}.${suffix_aligned}" : "${params.outdir}/preprocessing/splitncigar/${patient}/${id}/${cram.baseName}.${suffix_aligned}"
                crai = status <2 ? "${params.outdir}/preprocessing/mapped/${patient}/${id}/${crai.baseName.minus(".cram")}.${suffix_index}" : "${params.outdir}/preprocessing/splitncigar/${patient}/${id}/${crai.baseName.minus(".cram")}.${suffix_index}"
                table = "${params.outdir}/preprocessing/${patient}/${id}/recal_table/${id}.recal.table"
                ["sorted.csv", "ppatient,status,sample,lane,fastq_1,fastq_2,bam,bai,cram,crai,table,vcf\n${patient},${status},${sample},${lane},,,,,${cram},${crai},${table},\n"]
            }
        }
}
