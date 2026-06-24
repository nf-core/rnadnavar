process RUN_CONSENSUS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/00/0003a3962c416b4c458909d28e343758adf347e614507137b70ac5e4a3c094bb/data'
        : 'community.wave.seqera.io/library/bioconductor-complexheatmap_bioconductor-rtracklayer_r-data.table_r-dplyr_pruned:b2d4632cc9bf5c6e'}"

    input:
    tuple val(meta), path(input_file, stageAs: "inputs/*"), val(caller)

    output:
    tuple val(meta), path('*.consensus.vcf'), optional: true, emit: vcf
    tuple val(meta), path('*.consensus_*.vcf'), val(caller), optional: true, emit: vcf_separate
    tuple val(meta), path('*.consensus.maf'), optional: true, emit: maf
    tuple val(meta), path('*.consensus_*.maf'), val(caller), optional: true, emit: maf_separate
    path ("*.pdf"), optional: true, emit: pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    run_consensus.R --input_dir=inputs/ --out_prefix=${prefix}.consensus --cpu=${task.cpus} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version 2>&1) | head -n 1)
    END_VERSIONS
    """
}
