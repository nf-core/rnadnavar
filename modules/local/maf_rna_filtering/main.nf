process RNA_FILTERING {
    tag "$meta.id"
    label 'process_low'

    conda null
    container 'ghcr.io/raqmanzano/rnafilt:latest'

    input:
        tuple val(meta), path(maf), path(maf_realignment)
        tuple val(meta1), path(fasta)
        path fasta_fai

    output:
        tuple val(meta), path('*.maf'), emit: maf
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnadnavar/bin/
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def maf_realign_opt = maf_realignment? "--maf_realign $maf_realignment" : ""
        """
        filter_rna_mutations.py \\
            --maf $maf \\
            --ref $fasta \\
            --output ${prefix}.maf \\
            $maf_realign_opt \\
            $args
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(echo \$(python --version 2>&1) | sed 's/^.*Python (//;s/).*//')
        END_VERSIONS
        """

}
