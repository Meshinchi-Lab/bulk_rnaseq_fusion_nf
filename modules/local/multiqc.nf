process MULTIQC {

    //use image on quay.io
    container "quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0"

    // if process fails, retry running it
    errorStrategy "retry"

    input:
    path input
    val sample_sheet

    output:
    path("${sample_sheet}_multiqc_report.html")     , emit: report
    path("*report_data")                            , emit: data, optional: true

    script:
    def args = task.ext.args ?: ''
    """
    set -eou pipefail
    multiqc \\
        $args \\
        --filename "${sample_sheet}_multiqc_report.html" \\
        \$PWD
    """
}
