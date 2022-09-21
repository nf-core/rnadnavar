process DNA_FILTERING {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "anaconda::pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'quay.io/biocontainers/pandas:1.4.3' }"

    input:
        tuple val(meta), path(maf)
        path whitelist

    output:
        tuple val(meta), path('*.maf')                      , emit: maf
        path "versions.yml"                                 , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnadnavar/bin/
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        filter_dna_mutations.py -i $maf --output ${prefix}.maf $args
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python (//;s/).*//')
        END_VERSIONS
        """

}