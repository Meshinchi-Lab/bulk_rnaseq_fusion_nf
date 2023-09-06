process FASTQC {

    //use image on quay.io
    container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"

    input:
    tuple val(sample), file(R1), file(R2)

    output:
    path "${sample}", emit: fastqc

    script:
    def args = task.ext.args ?: ''
    """
    set -eou pipefail
    mkdir ${sample}
    fastqc \\
        --outdir ${sample} \\
        --threads ${task.cpus} \\
        --dir \$PWD \\
        $args \\
        $R1 $R2
    """
}
