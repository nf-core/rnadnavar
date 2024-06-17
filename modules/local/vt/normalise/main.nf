process VT_NORMALISE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::vt-0.57721-h17a1952_6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vt:0.57721--h17a1952_6' :
        'biocontainers/vt:0.57721--h17a1952_6' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta1), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "*.stats"                   , emit: stats
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vt \\
        normalize \\
        $vcf \\
        -r $fasta \\
        $args \\
        -o ${prefix}.vcf 2> ${prefix}.stats
    gzip ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vt normalize: \$(vt normalize -? 2>&1 | head -n1 | sed 's/^.*normalize //; s/ .*\$//')
    END_VERSIONS
    """












}
