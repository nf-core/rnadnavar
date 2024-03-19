//
// CHANNEL_CONSENSUS_CREATE_CSV
//

workflow CHANNEL_CONSENSUS_CREATE_CSV {
    take:
        maf_to_csv // channel: [mandatory] meta, maf, variantcaller
        csv_name

    main:
        // Creating csv files to restart from this step
        maf_to_csv.collectFile(keepHeader: true, skip: 1,sort: true, storeDir: "${params.outdir}/csv"){ meta, maf, variantcaller ->
            patient       = meta.patient
            sample        = meta.id
            status        = meta.status
            maf = "${params.outdir}/consensus/${variantcaller}/${meta.id}/${maf.getName()}"
            ["${csv_name}.csv", "patient,sample,status,variantcaller,maf\n${patient},${sample},${status},${variantcaller},${maf}\n"]
        }
}
