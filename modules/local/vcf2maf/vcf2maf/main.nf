process VCF2MAF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vcf2maf=1.6.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf2maf:1.6.21--hdfd78af_0' :
        'biocontainers/vcf2maf:1.6.21--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.maf")   , emit: maf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_decompressed = vcf.baseName.minus(".gz")
    def VERSION = '1.6.21' // Version information not provided by tool on CLI
    """
    gzip -d $vcf -c > ${vcf_decompressed}
    vcf2maf.pl \\
        --input-vcf ${vcf_decompressed} \\
        --output-maf ${prefix}.maf \\
        --ref-fasta $fasta \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION
    END_VERSIONS
    """
}
