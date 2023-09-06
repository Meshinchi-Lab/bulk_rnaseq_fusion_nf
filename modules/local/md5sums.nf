process MD5SUMS {

    container "centos/centos:centos7"

    // declare the input types and its variable names
    input:
    tuple val(meta), path(input)

    //define output files to save to the output_folder by publishDir command
    output:
    path "*.md5"

    script:
    def args = task.ext.args ?: ''
    """
    set -eou pipefail
    echo "Creating MD5sum checks"
    md5sum $args ${input} > ${input}.md5
    """
}
