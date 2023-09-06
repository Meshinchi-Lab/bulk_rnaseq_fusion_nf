nextflow.enable.dsl=2

//Define message for the process logs.
log.info """\
         R N A S E Q - F U S I O N  P I P E L I N E
         ===================================
         star_refs      :   ${params.star_genome_lib}
         cicero_refs    :   ${params.cicero_genome_lib}
         samples        :   ${params.sample_sheet}
         """
         .stripIndent()

include { star_index } from './subworkflows/local/star_index.nf'
include {
        MD5sums as md5_star
        MD5sums as md5_cicero } from './modules/local/fusion-processes.nf'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { 
        STAR_Prep_Fusion;
        STAR_Fusion; 
        fastqc; 
        multiqc; 
        STAR_aligner; 
        CICERO } from './modules/local/fusion-processes.nf'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:
nextflow run jennylsmith/STAR-fusion-NF <ARGUMENTS>
Required Arguments:
  Input Data:
  --sample_sheet        Single file with the location of all input data. Must be formatted
                        as a CSV with columns: Sample,R1,R2

  Required Genomic Reference Arguments:
  --star_genome_lib     The location of the CTAT Resource Library for STAR-Fusion - See https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/ for available ones
  --cicero_genome_lib   The location of the References for CICERO Fusion - See https://github.com/stjude/CICERO#reference for available references.
  --genome              The referenfence genome to be used by CICERO. A string value can be 'GRCh37-lite' or 'GRCh38_no_alt' only. 
  --build_index         Whether to build STAR aligner index. Boolean true or false. 

  Optional Arguments:
  --fasta_file          Location of directory which contains the reference genome fasta file. Required if build_index = true. 
  --gtf                 URL or path of the gtf file for the optional STAR index. Required if build_index = true. 
  
  Output Locations:
  --outdir              The location of the pipeline results files 
  --star_index_dir      The location of the STAR index for running STAR aligner 
 """.stripIndent()
}

workflow  fusion_calls {
    // Define the input paired fastq files in a sample sheet and genome references.
    // The sample_sheet is comma separated with column names "Sample","R1","R2"
    Channel.fromPath(file(params.sample_sheet))
        .splitCsv(header: true, sep: ',')
        .map { sample -> [ sample["Sample"] + "_", 
                           file(sample["R1"], checkIfExists: true), 
                           file(sample["R2"], checkIfExists: true) ]
            }
        .set { fqs_ch }

    // QC on the fastq files
    fastqc(fqs_ch)

    // Prepare chimeric junctions files and input into STAR-fusion
    Channel.fromPath(file(params.star_genome_lib, checkIfExists: true))
        .collect()
        .set { star_genome_lib }
    STAR_Prep_Fusion(star_genome_lib, fqs_ch)
    STAR_Fusion(star_genome_lib, fqs_ch, STAR_Prep_Fusion.out.chimera)
    md5_star(STAR_Prep_Fusion.out.bam)

    // CICERO requires GRCh37-lite or GRCh38_no_alt aligned BAMs,
    // STAR-aligner must be re-run on these specific assemblies 
    if ( params.build_index ) {
        // optionally, build the index from fasta and gtf
        star_index()
        star_index.out.index
            .collect()
            .set { star_index }
    } else {
        // Stage the index files
        Channel.fromPath(file(params.star_index_out, checkIfExists: true))
            .collect()
            .set { star_index }
    }
    star_index.view {" the index channel is $it"}
    STAR_aligner(star_index, fqs_ch)
    SAMTOOLS_INDEX(STAR_aligner.out.bam)
    md5_cicero(STAR_aligner.out.bam)
    STAR_aligner.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set { bam_bai_ch }

    // Run CICERO on the STAR-aligned BAM files.
    Channel.fromPath(file(params.cicero_genome_lib, checkIfExists: true))
        .collect()
        .set { cicero_genome_lib }
    Channel.value(params.genome)
        .set { genome }
    CICERO(bam_bai_ch, cicero_genome_lib, genome)

    // MultiQC 
    def sample_sheet = file(params.sample_sheet).simpleName
    fastqc.out.fastqc
        .collect()
        .set { mqc_ch }
    multiqc(mqc_ch, sample_sheet)
}