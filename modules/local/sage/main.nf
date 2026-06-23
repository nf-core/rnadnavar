process SAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hmftools-sage=3.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.2.3--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-sage:3.2.3--hdfd78af_0' }"

    input:
    tuple val(meta),  path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor), path(intervals)
    tuple val(meta1), path(ensembl_dir)
    tuple val(meta2), path(sage_high_confidence)
    tuple val(meta3), path(sage_actionable_panel)
    tuple val(meta4), path(sage_known_hotspots)
    tuple val(meta5), path(fasta)
    tuple val(meta6), path(fai)
    tuple val(meta7), path(dict)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args   ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def reference        = input_normal    ? "-reference ${meta.normal_id} -reference_bam ${input_normal}" : ""
    def interval_prep    = intervals       ? "INTER=\$(sed -E 's/\\s+0\\s+/\\t1\\t/g' $intervals | sed 's/\t/:/g' | paste -s -d ';')" : ""
    def interval_arg     = intervals       ? "-specific_regions \$INTER" : ""
    def interval_notice  = intervals       ? "echo \"[WARNING] If no reads in the intervals from $intervals, sage won't work\"" : ""
    def ensembl_data_dir = ensembl_dir     ? "-ensembl_data_dir ${ensembl_dir}": ""
    def avail_mem        = 3072
    if (!task.memory) {
        log.info '[SAGE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    ${interval_notice}
    export _JAVA_OPTIONS="-Xmx${avail_mem}g"
    ${interval_prep}

    SAGE \\
        -out ${prefix}.vcf \\
        -ref_genome $fasta \\
        -threads $task.cpus \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${input_tumor} \\
        -high_confidence_bed ${sage_high_confidence} \\
        -panel_bed ${sage_actionable_panel} \\
        -hotspots ${sage_known_hotspots} \\
        $interval_arg \\
        $reference \\
        $ensembl_data_dir \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SAGE: \$(SAGE 2>&1 | grep 'Sage version' | cut -d " " -f6)
    END_VERSIONS
    """
    stub:

    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SAGE: \$(SAGE 2>&1 | grep 'Sage version' | cut -d " " -f6)
    END_VERSIONS
    """

}
