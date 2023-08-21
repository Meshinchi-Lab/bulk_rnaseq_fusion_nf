process MULTIQC {

    publishDir "$params.multiQC"

    //use image on quay.io
    container "quay.io/lifebitai/MULTIQC:latest"
    cpus 2
    memory "16 GB"

    // if process fails, retry running it
    errorStrategy "retry"

    input:
    path "*"
    val sample_sheet

    output:
    path "${sample_sheet}_MULTIQC_report.html"

    script:
    """
    set -eou pipefail

    MULTIQC -v --filename "${sample_sheet}_multiqc_report.html"  .
    """
}
