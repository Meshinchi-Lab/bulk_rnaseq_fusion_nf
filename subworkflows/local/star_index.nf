include {
        UNZIP as gunzip_fasta
        UNZIP as gunzip_gtf } from '../../modules/local/unzip.nf'
include { STAR_INDEX } from '../../modules/local/star/star_index.nf'

//Workflow to create an index for STAR-aligner given a human genome fasta and the location of the GTF file. 
workflow star_index {

    main:
    //Stage the genome fasta files for the index building step
    Channel.fromPath(file(params.fasta_file, checkIfExists: true))
        .set{ fasta }
    // if fasta  is gzipped, it must be decompressed   
    if(params.fasta_file.endsWith(".gz")){
        gunzip_fasta(fasta)
        gunzip_fasta.out.file
            .set{ fasta }
    }
    //Download and stage the GTF file 
    Channel.fromPath(file(params.gtf, checkIfExists: true))
        .set{ gtf }
     //if gtf is gzipped, it must be decompressed   
    if(params.gtf.endsWith(".gz")){
        gunzip_gtf(gtf)
        gunzip_gtf.out.file
            .set{ gtf }
    } 
    //STAR genomeGenerate to build the index
    STAR_INDEX(fasta, gtf)

    emit:
    index           = STAR_INDEX.out.index
    fasta           = fasta
    gtf             = gtf
}