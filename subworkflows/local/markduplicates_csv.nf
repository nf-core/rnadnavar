//
// MARKDUPLICATES_CSV
//

workflow MARKDUPLICATES_CSV {
    take:
        cram_markduplicates // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_markduplicates.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, file, index ->
            id             = meta.id
            lane           = meta.lane
            patient        = meta.patient
            sample         = meta.sample
            status         = meta.status
            suffix_aligned = params.save_output_as_bam ? "bam" : "cram"
            suffix_index   = params.save_output_as_bam ? "bam.bai" : "cram.crai"
            file   = "${params.outdir}/preprocessing/markduplicates/${patient}/${id}/${file.baseName}.${suffix_aligned}"
            index   = "${params.outdir}/preprocessing/markduplicates/${patient}/${id}/${index.baseName.minus(".cram")}.${suffix_index}"
            ["markduplicates_no_table.csv", "patient,status,sample,lane,fastq_1,fastq_2,bam,bai,cram,crai,table,vcf\n${patient},${status},${sample},${lane},,,,,${file},${index},,\n"]
        }
}
