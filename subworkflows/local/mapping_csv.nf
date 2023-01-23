//
// MAPPING_CSV
//

workflow MAPPING_CSV {
    take:
        bam_indexed         // channel: [mandatory] meta, bam, bai

    main:
        // Creating csv files to restart from this step
        bam_indexed.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, bam, bai ->
            id = meta.id
            patient = meta.patient
            sample  = meta.sample
            status  = meta.status
            lane    = meta.lane
            bam     = "${params.outdir}/preprocessing/mapped/${patient}/${id}/${bam.name}"
            bai     = "${params.outdir}/preprocessing/mapped/${patient}/${id}/${bai.name}"
            ["mapped.csv", "patient,status,sample,lane,fastq_1,fastq_2,bam,bai,cram,crai,table,vcf\n${patient},${status},${sample},${lane},,,${bam},${bai},,,,\n"]
        }
}
