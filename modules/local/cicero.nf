process CICERO {

    // use CICERO container from stjude
    container "ghcr.io/stjude/cicero:v1.9.6"

    // declare the input types and its variable names
    input:
    tuple val(sample), file(BAM), file(BAI)
    path cicero_genome_lib
    val genome

    //define output files
    output:
    path("${sample}/*/*final_fusions.txt")  , emit: cicero, optional: true
    path("${sample}/*/*.txt")               , emit: outfiles
    path("${sample}/*.log")                 , emit: logs

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -eou pipefail
    export TMPDIR="\$PWD"

    # CICERO fusion detection algorithm
    Cicero.sh \\
        -n ${task.cpus} \\
        -b \$PWD/$BAM \\
        -g ${genome} \\
        -r \$PWD/$cicero_genome_lib \\
        $args \\
        -o ${sample}
    """
}
