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
    path("${sample}/CICERO_DATADIR/${BAM.baseName}/*.txt")  , emit: outfiles
    path("${sample}/*.{err,log,out}")                       , emit: logs

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample}"
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TMPDIR=\$PWD

    Cicero.sh \\
        -n ${task.cpus} \\
        -b \${PWD}/$BAM \\
        -g ${genome} \\
        -r \${PWD}/$cicero_genome_lib \\
        $args \\
        -o ${sample}
    """
}
