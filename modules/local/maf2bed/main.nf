process MAF2BED {
    tag "$meta.id"
    label 'process_single'

    conda "anaconda::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    tuple val(meta), path(maf)

    output:
    tuple val(meta), path('*.bed') , emit: bed
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
#!/usr/bin/env python
import subprocess
import pandas as pd
maf = pd.read_csv("${maf}", sep="\\t", comment="#")
bed = maf[["Chromosome", "Start_Position", "End_Position"]]
bed.to_csv("${prefix}.bed", sep="\\t", index=False, header=False)
subprocess.check_call('cat <<-END_VERSIONS > versions.yml\\n\\"${task.process}\\":\\n\tpython: \$(echo \$(python --version 2>&1) | sed \\"s/^.*Python (//;s/).*//\\")\\nEND_VERSIONS', shell=True)
    """
}
