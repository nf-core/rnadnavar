//
// RECALIBRATE_CSV
//

workflow RECALIBRATE_CSV {
    take:
        cram_recalibrated_index // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_recalibrated_index.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, file, index ->
            patient = meta.patient
            sample  = meta.sample
            status  = meta.status
            id      = meta.id
            lane    = meta.lane
            file = "${params.outdir}/preprocessing/recalibrated/${patient}/${id}/${file.name}"
            index = "${params.outdir}/preprocessing/recalibrated/${patient}/${id}/${index.name}"
            ["recalibrated.csv","patient,status,sample,lane,fastq_1,fastq_2,bam,bai,cram,crai,table,vcf\n${patient},${status},${sample},${lane},,,,,${file},${index},,\n"]
        }
}
