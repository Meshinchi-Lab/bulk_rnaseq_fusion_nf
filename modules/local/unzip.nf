process UNZIP {

    container "centos/centos:centos7"

    input:
    path zipped_file  

    output:
    path("*"), emit: file

    script:
    def args = task.ext.args ?: ''
    def file_ext = "${zipped_file.getClass()}"
    // def outfile = "${zipped_file.take())}"
    """
    echo "the outfile is $file_ext"
    gunzip $args -f $zipped_file 
    """
}
