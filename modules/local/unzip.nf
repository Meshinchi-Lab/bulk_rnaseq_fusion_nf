process UNZIP {

    container "centos/centos:centos7"

    input:
    path zipped_file  

    output:
    path("*"), emit: file

    script:
    def args = task.ext.args ?: ''
    def outfile = "${zipped_file.toString()}"
    // def outfile = "${zipped_file.take())}"
    """
    echo "the outfile is $outfile"
    gunzip $args -f $zipped_file 
    """
}
