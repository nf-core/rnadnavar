process RNA_FILTERING {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "anaconda::pandas=1.4.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2%3A3b083bb5eae6e491b8579589b070fa29afbea2a1-0' :
        'quay.io/biocontainers/mulled-v2-629aec3ba267b06a1efc3ec454c0f09e134f6ee2%3A3b083bb5eae6e491b8579589b070fa29afbea2a1-0' }"
    input:
        tuple val(meta), path(maf_first_pass), path(maf_second_pass)
        path fasta

    output:
        tuple val(meta), path('*.maf')                     , emit: maf
        path "versions.yml"                                , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnadnavar/bin/
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def maf_second = maf_second_pass.name != 'NO_FILE.maf' ? "--maf_2pass $maf_second_pass" : ''
        """
        filter_rna_mutations.py \\
            --maf $maf_first_pass \\
            --ref $fasta \\
            --output ${prefix}.maf \\
            $maf_second \\
            $args
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python (//;s/).*//')
        END_VERSIONS
        """

}