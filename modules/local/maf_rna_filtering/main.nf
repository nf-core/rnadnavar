process RNA_FILTERING {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/raqmanzano/rnafilt:latest' : null }"

    input:
        tuple val(meta), path(maf_first_pass), path(maf_second_pass)
        path fasta
        path fasta_fai

    output:
        tuple val(meta), path('*.maf')                     , emit: maf
        path "versions.yml"                                , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnadnavar/bin/
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def maf_second_opt = maf_second_pass? "--maf_2pass $maf_second_pass" : ""
        """
        filter_rna_mutations.py \\
            --maf $maf_first_pass \\
            --ref $fasta \\
            --output ${prefix}.maf \\
            $maf_second_opt \\
            $args
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python (//;s/).*//')
        END_VERSIONS
        """

}
