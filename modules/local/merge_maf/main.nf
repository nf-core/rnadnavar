process MERGE_MAF {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(input_files, stageAs: "inputs/*")

    output:
    tuple val(meta), path("*merged.maf"), emit: maf_out
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Get header(s)
    ls inputs/* |head -n1 | xargs grep -m 2 -E 'Hugo_Symbol|^#' > ${prefix}.merged.maf
    # Concat files without header
    grep -hEv 'Hugo_Symbol|^#' inputs/* >> ${prefix}.merged.maf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$( cat --version 2>&1 | head -n 1 )
    END_VERSIONS
    """
}
