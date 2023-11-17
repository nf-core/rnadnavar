//
// CHANNEL_VARIANT_CALLING_CREATE_CSV
//

workflow CHANNEL_VARIANT_CALLING_CREATE_CSV {
    take:
        vcf_to_csv // channel: [mandatory] meta, vcf
        csv_name

    main:
        // Creating csv files to restart from this step
        vcf_to_csv.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, vcf ->
            patient       = meta.patient
            sample        = meta.id
            variantcaller = meta.variantcaller
            status        = meta.status
            vcf = "${params.outdir}/variant_calling/${variantcaller}/${meta.id}/${vcf.getName()}"
            ["${csv_name}.csv", "patient,sample,status,variantcaller,vcf\n${patient},${sample},${status},${variantcaller},${vcf}\n"]
        }
}
