include {
        UNZIP as gunzip_fasta;
        UNZIP as gunzip_gtf } from '../../modules/local/unzip.nf'
include { BUILD_GENOME_REFS } from '../../modules/local/build_genome_refs.nf'

//Workflow to create an index and fusion genome library for STAR-aligner given a genome fasta and the location of the GTF file. 
workflow build_genome_lib {
    take:
    fasta
    gtf

    main:
    if ( fasta instanceof String ){
        //Stage the genome fasta files for the index building step
        Channel.fromPath(file(fasta, checkIfExists: true))
                .set{ fasta }
        // if fasta  is gzipped, it must be decompressed   
        if(params.fasta_file.endsWith(".gz")){
            gunzip_fasta(fasta)
            gunzip_fasta.out.file
                .set{ fasta }
        }
    }

    // Stage the GTF file 
    if ( gtf instanceof String ){
        Channel.fromPath(file(gtf, checkIfExists: true))
            .set{ gtf }
        //if gtf is gzipped, it must be decompressed   
        if(params.gtf.endsWith(".gz")){
            gunzip_gtf(gtf)
            gunzip_gtf.out.file
                .set{ gtf }
        }
    }

    // prepare pfam database from path or allow prep_genome_lib.pl to download the current DB
    if ( params.pfam_db.matches('current') ){
        Channel.value(params.pfam_db)
            .set { pfam }
    } else {
        Channel.fromPath(file(params.pfam_db, checkIfExists: true))
            .set { pfam }
    }

    // prepare dfam database from path or allow prep_genome_lib.pl to download the human/mouse DB
    if ( params.dfam_db.matches("human|mouse") ){
         Channel.value(params.dfam_db)
            .set { dfam }
    } else {
        Channel.fromPath(file(params.dfam_db, checkIfExists: true))
            .set { pfam }
    }

    // STAR prep_genome_lib.pl to build the index and fusion refs
    Channel.fromList( [ fasta, gtf, dfam, pfam ] )
        .collect()
        .set { ref_files }
    ref_files.view {"the ref files are $it"}
    BUILD_GENOME_REFS(ref_files)

    emit:
    genome_lib           = BUILD_GENOME_REFS.out.genome_dir
}