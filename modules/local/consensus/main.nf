process RUN_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-rtracklayer bioconda::bioconductor-complexheatmap conda-forge::r-ggrepel conda-forge::r-data.table conda-forge::r-dplyr conda-forge::ggpubr "
    container 'ghcr.io/raqmanzano/renv:latest'

    input:
        tuple val(meta), path(vcf, stageAs: "inputs/*"), val(caller)

    output:
        tuple val(meta), path('*.consensus.vcf')                , optional:true , emit: vcf
        tuple val(meta), path('*.consensus_*.vcf'), val(caller) , optional:true , emit: vcf_separate
        tuple val(meta), path('*.consensus.maf')                , optional:true , emit: maf
        tuple val(meta), path('*.consensus_*.maf'), val(caller) , optional:true , emit: maf_separate
        path("*.pdf")                                           , optional:true , emit: pdf
        path "versions.yml"                                                     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/rnadnavar/bin/
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        run_consensus.R --input_dir=inputs/ --out_prefix=${prefix}.consensus --cpu=$task.cpus $args
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$(echo \$(R --version 2>&1) | head -n 1)
        END_VERSIONS
        """
}
