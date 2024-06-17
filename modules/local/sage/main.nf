process SAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hmftools-sage=3.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.2.3--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-sage:hmftools-sage:3.2.3--hdfd78af_0' }"

    input:
    tuple val(meta),  path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor), path(intervals)
    tuple val(meta1), path(ensembl_dir)
    tuple val(meta2), path(sage_highconfidence)
    tuple val(meta3), path(sage_actionablepanel)
    tuple val(meta4), path(sage_knownhotspots)
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
    def interval_cmd     = intervals       ? "INTER=\$(sed -E 's/\\s+0\\s+/\\t1\\t/g' $intervals | sed 's/\t/:/g' | paste -s -d ';')": ""
    def interval         = intervals       ? "-specific_regions \$INTER": ""
    def ensembl_data_dir = ensembl_dir     ? "-ensembl_data_dir ${ensembl_dir}": ""
    def avail_mem        = 3072
    if (!task.memory) {
        log.info '[SAGE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    if (intervals){
    """
    echo "[WARNING] If no reads in the intervals from $intervals, sage won't work"
    export _JAVA_OPTIONS="-Xmx${avail_mem}g"
    ${interval_cmd}

    SAGE \\
        -out ${prefix}.vcf \\
        -ref_genome $fasta \\
        -threads $task.cpus \\
        -tumor ${meta.tumor_id} -tumor_bam ${input_tumor} \\
        -high_confidence_bed ${sage_highconfidence} \\
        -panel_bed ${sage_actionablepanel} \\
        -hotspots ${sage_knownhotspots} \\
        $reference \\
        $interval \\
        $ensembl_data_dir \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SAGE: \$(grep 'Sage version' .command.log | cut -d " " -f6)
    END_VERSIONS
    """

    } else {
    """
    export _JAVA_OPTIONS="-Xmx${avail_mem}g"
    SAGE \\
        -out ${prefix}.vcf \\
        -ref_genome $fasta \\
        -threads $task.cpus \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${input_tumor} \\
        -high_confidence_bed ${sage_highconfidence} \\
        -panel_bed ${sage_actionablepanel} \\
        -hotspots ${sage_knownhotspots} \\
        $ensembl_data_dir \\
        $reference \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SAGE: \$(grep 'Sage version' .command.log | cut -d " " -f6)
    END_VERSIONS
    """
    }
    stub:

    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SAGE: \$(SAGE | grep 'Sage version' .command.log | cut -d " " -f6)
    END_VERSIONS
    """

}
