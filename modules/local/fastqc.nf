process FASTQC {

    //use image on quay.io
    container "quay.io/biocontainers/FASTQC:0.11.9--hdfd78af_1"
    cpus 2
    memory "16 GB"

    // if process fails, retry running it
    errorStrategy "retry"

    input:
    tuple val(Sample), file(R1), file(R2)

    output:
    path "FASTQC_${Sample}_logs"

    script:
    """
    mkdir FASTQC_${Sample}_logs
    FASTQC -o fastqc_${Sample}_logs -t 6 -f fastq -q $R1 $R2
    #rm $R1 $R2 to avoid repetitive upload to the S3 bucket
    """
}
