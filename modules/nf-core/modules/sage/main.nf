def VERSION = '3.1'  // Version information not provided by tool on CLI

process SAGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::hmftools-sage=3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:3.1--hdfd78af_0' :
        'quay.io/biocontainers/hmftools-sage:3.1--hdfd78af_0' }"

    input:
        tuple val(meta), path(normal), path(normal_index), path(tumor), path(tumor_index), path(intervals)
        path fasta
        path fasta_fai
        path dict
        path highconfidence
        path actionablepanel
        path knownhot
        path ensbl_sage

    output:
        tuple val(meta), path("*.vcf"), emit: vcf
        path  "versions.yml"          , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        def reference             = normal          ? "-reference ${meta.normal_id} -reference_bam ${normal}" : ""
        def HighConfidence        = highconfidence  ? "-high_confidence_bed ${highconfidence}"                : ""
        def ActionableCodingPanel = actionablepanel ? "-panel_bed ${actionablepanel}"                         : ""
        def KnownHotspots         = knownhot        ? "-hotspots ${knownhot}"                                 : ""
        def avail_mem = 4
        if (!task.memory) {
            log.info '[SAGE] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = task.memory.giga
        }
        if (intervals){ // If no reads the intervals don't work in sage
        """
        export _JAVA_OPTIONS="-Xmx${avail_mem}g"
        INTER=\$(sed -E 's/\\s+0\\s+/\\t1\\t/g' $intervals | grep -v chrM | sed 's/\t/:/g' | paste -s -d ';')

        SAGE \\
            -out ${prefix}.vcf \\
            -ref_genome $fasta \\
            -threads $task.cpus \\
            -tumor ${meta.tumor_id} \\
            -tumor_bam ${tumor} \\
            $reference \\
            -ensembl_data_dir $ensbl_sage \\
            $HighConfidence \\
            $ActionableCodingPanel \\
            $KnownHotspots \\
            -specific_regions \$INTER \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sage: $VERSION
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
            -tumor_bam ${tumor} \\
            $reference \\
            -ensembl_data_dir $ensbl_sage \\
            $HighConfidence \\
            $ActionableCodingPanel \\
            $KnownHotspots \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sage: $VERSION
        END_VERSIONS
        """
        }



}
