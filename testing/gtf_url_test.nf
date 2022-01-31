nextflow.enable.dsl=2
params.gtf_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
params.genome_fasta="data"
// params.genome_fasta="gzipped_data"

//Must include the same module with different names to be used in the same workflow
include {
    unzip as gunzip_fasta
    unzip as gunzip_gtf } from './modules/test_modules.nf'

include { echo_file } from './modules/test_modules.nf'

//https://stackoverflow.com/questions/66654223/nextflow-groovy-check-if-item-in-channel-list
//https://github.com/nextflow-io/patterns/blob/master/docs/process-when-empty.adoc
// check_gz = fasta.filter( ~/^.+gz/ ).ifEmpty{ 'EMPTY' }
def check(x) {
    result = x ==~ 'EMPTY'
    return(result)
}

workflow  {
  //Download and stage the GTF file from a given URL 
  Channel.fromPath(params.gtf_url)
      .ifEmpty { error  "No file found at URL ${params.gtf_url}" }
      .set{gtf}       
  if(params.gtf_url.endsWith(".gz")){
    // echo_file(file(params.gtf_url))
    //   .set{gtf}
    gunzip_gtf(gtf)
    gunzip_gtf.out.unzipped_file.set{gtf}
  } 
  gtf.view()

  //Stage the genome fasta files for the index building step
  fastq_pattern = "${params.genome_fasta}/*.{fa,fasta,fa.gz,fasta.gz}"
  Channel.fromPath(fastq_pattern)
    .ifEmpty { error "No files found matching the pattern ${fastq_pattern}" }
    .set{fasta}

  //GOAL: if channel contains gzipped files, then unzip them. 
  check = fasta.filter( ~/^.+gz/ ).ifEmpty{ 'EMPTY' } =~ 'EMPTY' ? "not_gzipped" : "gzipped"
  // println fasta.filter( ~/^.+gz/ ).ifEmpty{ 'EMPTY' }.toString() ==~ 'EMPTY' //FALSE
  // println 'EMPTY' ==~ 'EMPTY' //TRUE 
  //  println "the files are: ${check}"  

  x = fasta.filter( ~/^.+gz/ ).ifEmpty{ 'EMPTY' }
  x.view()
  y = check(x)
  println "the output is ${y}"

  // x = 'foo'
  // if( x ==~ /foo/ ){
  //   println "foo works"
  // }

  // if( check ==~ /gzipped/ ){
  //   gunzip_fasta(fasta)
  //   gunzip_fasta.out.unzipped_file.set{fasta}
  // }
  // fasta.view()
}

