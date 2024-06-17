process VCFFILTER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcftools=0.1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcflib:1.0.3--hecb563c_1' :
        'biocontainers/vcflib:1.0.3--hecb563c_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    gunzip -c $vcf | vcffilter \\
        $args > ${prefix}.vcf
    gzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcffilter: \$(echo \$(vcffilter 2>&1) | sed 's/^.*vcflib ( //;s/).*//' | cut -f2 -d ' ')
    END_VERSIONS
    """
}
