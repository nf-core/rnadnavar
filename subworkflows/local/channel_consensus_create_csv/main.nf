//
// CHANNEL_CONSENSUS_CREATE_CSV
//

workflow CHANNEL_CONSENSUS_CREATE_CSV {
    take:
        maf_to_csv // channel: [mandatory] [meta, maf] or [meta, maf, variantcaller]
        csv_name

    main:
        // Creating csv files to restart from this step
        // Consensus/rescue outputs now flow downstream as canonical [meta, maf] tuples.
        // For restart CSVs we still need a variantcaller column, so infer it here when absent.
        maf_to_csv.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/csv") { row ->
            def meta = row[0]
            def maf = row[1]
            def variantcaller = meta.variantcaller
            def patient       = meta.patient
            def sample        = meta.id
            def status        = meta.status
            def subfolder = maf.getName().contains('consensus') ? 'consensus' : 'vcf2maf'
            maf = "${params.outdir}/consensus/${subfolder}/${meta.id}/${maf.getName()}"
            ["${csv_name}.csv", "patient,sample,status,variantcaller,maf\n${patient},${sample},${status},${variantcaller},${maf}\n"]
        }
}
