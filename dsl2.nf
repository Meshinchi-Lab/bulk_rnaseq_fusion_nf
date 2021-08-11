params.reads = "$baseDir/data/ggal/gut_{1,2}.fq"

Channel
    .from( 1, 2, 3, 4 )
    .collect()
    .view()

Channel.fromFilePairs( params.reads, checkIfExists: true ).view()
