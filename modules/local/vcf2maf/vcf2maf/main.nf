process VCF2MAF {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vcf2maf-1.6.21-hdfd78af_0" : null )
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vcf2maf:1.6.21--hdfd78af_0' :
        'quay.io/biocontainers/vcf2maf:1.6.21--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)
    path(fasta)
    val(genome)

    output:
    tuple val(meta), path("*.maf.gz"), emit: maf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_decompressed = vcf.baseName.minus(".gz")
    """
    gzip -d $vcf -c > ${vcf_decompressed}.vcf
    vcf2maf.pl \\
        --input-vcf ${vcf_decompressed}.vcf \\
        --output-maf ${prefix}.maf \\
        --inhibit-vep \\
        --ref-fasta $fasta \\
        --normal-id $meta.normal_id \\
        --tumor-id $meta.tumor_id \\
        --vcf-tumor-id $meta.tumor_id \\
        --vcf-normal-id $meta.normal_id \\
        --max-subpop-af 0.0001 \\
        --retain-fmt AD,DP,AF \\
        $args



    gzip ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: \$(vcf2maf -v? 2>&1 | head -n1 | sed 's/^.*vcf2maf //; s/ .*\$//')
    END_VERSIONS
    """

}