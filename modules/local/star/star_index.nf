process STAR_INDEX {

    // use image on quay.io
    // container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"

    //input genome fasta and gtf
    input: 
    path fasta
    path gtf
    
    //output the index into a diretory, and the logfile
    output:
    path("GenomeDir"), emit: index

    script:
    def args = task.ext.args ?: ''
    """
    mkdir \$PWD/GenomeDir
    STAR --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        $args \\
        --genomeDir \$PWD/GenomeDir \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf
    """
}
