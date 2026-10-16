//
// ANNOTATION CACHE INITIALISATION
//

// Initialise channels based on params or indices that were just built
// For all modules here:
// A when clause condition is defined in the conf/modules.config to determine if the module should be run
// Condition is based on params.step and params.tools
// If and extra condition exists, it's specified in comments

workflow ANNOTATION_CACHE_INITIALISATION {
    take:
    vep_enabled
    vep_cache
    vep_species
    vep_cache_version
    vep_genome
    vep_custom_args
    help_message

    main:

    if (vep_enabled) {
        def vep_annotation_cache_key = (vep_cache == "s3://annotation-cache/vep_cache/") ? "${vep_cache_version}_${vep_genome}/" : ""
        def vep_species_suffix = vep_custom_args.contains("--merged") ? '_merged' : (vep_custom_args.contains("--refseq") ? '_refseq' : '')
        def vep_cache_dir = "${vep_annotation_cache_key}${vep_species}${vep_species_suffix}/${vep_cache_version}_${vep_genome}"
        def vep_cache_path_full = file("$vep_cache/$vep_cache_dir", type: 'dir')
        if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
            if (vep_cache == "s3://annotation-cache/vep_cache/") {
                error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
            } else {
                error("Path provided with VEP cache is invalid.\nMake sure there is a directory named ${vep_cache_dir} in ${vep_cache}./n${help_message}")
            }
        }
        // Keep the cache channel aligned with the nf-core VEP module contract:
        // a single reusable value channel carrying [meta, cache_dir].
        ensemblvep_cache = Channel.value([[id: vep_cache_path_full.baseName], vep_cache_path_full])
    } else ensemblvep_cache = []

    emit:
    ensemblvep_cache // channel: [ meta, cache ]
}
