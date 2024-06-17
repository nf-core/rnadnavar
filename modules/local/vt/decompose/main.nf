process VT_DECOMPOSE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::vt-0.57721-h17a1952_6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vt:0.57721--h17a1952_6' :
        'biocontainers/vt:0.57721--h17a1952_6' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "*.stats"                   , emit: stats
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_decompressed = vcf.baseName.minus(".gz")
    """
    gzip -d $vcf -c > ${vcf_decompressed}
    # GQ is a float when empty which can happen with some tools like freebayes - this is a fix
    sed -i -E 's/(##FORMAT=<ID=GQ\\S+)(Integer)/\\1Float/' ${vcf_decompressed}

        vt \\
            decompose \\
            ${vcf_decompressed} \\
            $args \\
            -o ${prefix}.vcf 2> ${prefix}.stats
        gzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt decompose: \$(vt decompose -? 2>&1 | head -n1 | sed 's/^.*decompose //; s/ .*\$//')
    END_VERSIONS
    """
}
