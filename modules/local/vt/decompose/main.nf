process VT_DECOMPOSE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "" : null )
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vt:0.57721--h17a1952_6' :
        'quay.io/biocontainers/vt:0.57721--h17a1952_6' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vt \\
        decompose \\
        $vcf \\
        $args \\
        -o ${prefix}.vcf 2> ${prefix}.stats
    gzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt decompose: \$(vt decompose -? 2>&1 | head -n1 | sed 's/^.*decompose //; s/ .*\$//')
    END_VERSIONS
    """












}