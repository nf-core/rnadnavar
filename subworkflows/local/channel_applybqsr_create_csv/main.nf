//
// CHANNEL_APPLYBQSR_CREATE_CSV
//

workflow CHANNEL_APPLYBQSR_CREATE_CSV {
    take:
        cram_recalibrated_index // channel: [mandatory] meta, cram, crai

    main:
        // Creating csv files to restart from this step
        cram_recalibrated_index.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { meta, file, index ->
            def patient = meta.patient
            def sample  = meta.sample
            def status  = meta.status
            file = "${params.outdir}/preprocessing/recalibrated/${sample}/${file.name}"
            index = "${params.outdir}/preprocessing/recalibrated/${sample}/${index.name}"

            def type = params.save_output_as_bam ? "bam" : "cram"
            def type_index = params.save_output_as_bam ? "bai" : "crai"

            ["recalibrated.csv", "patient,status,sample,${type},${type_index}\n${patient},${status},${sample},${file},${index}\n"]
        }
}
