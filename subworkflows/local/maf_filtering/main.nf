//
// Filtering somatic mutation analysis
//
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
include { MAF_FILTERING as FILTERING } from '../../../modules/local/maf_filtering/main'

workflow MAF_FILTERING {

    take:
    maf_to_filter
    fasta
    input_sample

    main:
    versions  = Channel.empty()
    maf       = Channel.empty()
    if (params.step in ['mapping', 'markduplicates', 'splitncigar',
                'prepare_recalibration', 'recalibrate', 'variant_calling',
                'normalise', 'consensus', 'filtering'] ) {

        if (params.step == 'filtering') maf_to_filter = input_sample
		maf_to_filter.dump(tag:"maf_to_filter")
        // BASIC FILTERING
        FILTERING(maf_to_filter, fasta)
		maf      = FILTERING.out.maf
        versions = versions.mix(FILTERING.out.versions)
	}

    emit:
    maf        = maf
    versions   = versions                                                         // channel: [ versions.yml ]
}