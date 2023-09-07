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

// Subworkflows
include { star_index } from './subworkflows/local/star_index.nf'
include { build_genome_lib } from './subworkflows/local/star_build_refs.nf'

// QC modules
include {
        MD5SUMS as md5_star;
        MD5SUMS as md5_cicero } from './modules/local/md5sums.nf'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { FASTQC } from './modules/local/fastqc.nf'
include { MULTIQC } from './modules/local/multiqc.nf'

// Fusion modules
include { STAR_FUSION } from './modules/local/star_fusion.nf'
include { STAR_PREP_FUSION } from './modules/local/star_prep_fusion.nf'
include { STAR_ALIGNER } from './modules/local/star_aligner.nf'
include { CICERO } from './modules/local/cicero.nf'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:
nextflow run RSC-RP/star_fusion_nf <ARGUMENTS>
  Required Arguments:
  --sample_sheet        Single file with the location of all input data. Must be formatted as a CSV with columns: Sample,R1,R2
  --project             Cybertron project code (eg 207f23bf-acb6-4835-8bfe-142436acb58c) for HPC.
  --queue               Cybertron PBS queue name (eg paidq) for HPC.

  Required Genomic Reference Arguments:
  --star_genome_lib     The location of the CTAT Resource Library for STAR-Fusion - See https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/ for available ones
  --cicero_genome_lib   The location of the References for CICERO Fusion - See https://github.com/stjude/CICERO#reference for available references.
  --genome              The referenfence genome to be used by CICERO. A string value can be 'GRCh37-lite' or 'GRCh38_no_alt' only. 
  --build_index         Whether to build STAR aligner index. Boolean true or false. 
  --star_index_dir      The location of the STAR index for running STAR aligner. Can exist already or it will be created when build_index = true.

  Optional Arguments:
  --fasta_file          Location of directory which contains the reference genome fasta file. Required if build_index = true. 
  --gtf                 URL or path of the gtf file for the optional STAR index. Required if build_index = true. 
  
  Output Locations:
  --outdir              The location of the pipeline results files.
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
    FASTQC(fqs_ch)

    // Prepare STAR aligner index
    def build_index = params.build_index.toBoolean()
    def build_genome_lib = params.build_genome_library.toBoolean()
    if ( build_index ) {
        // optionally, build the index from fasta and gtf
        star_index()
        star_index.out.index
            .collect()
            .set { star_index }
    } else {
        // Stage the index files
        Channel.fromPath(file(params.star_index_dir, checkIfExists: true))
            .collect()
            .set { star_index }
    }
    // Optionally, Prepare STAR-fusion references
    if ( build_genome_lib ){
        if ( build_index ) {
            build_genome_lib(star_index.out.fasta, star_index.out.gtf)
        } else {
            build_genome_lib(params.fasta_file, params.gtf)
        }
        build_genome_lib().out.genome_lib
            .set { star_genome_lib }
    } else {
        Channel.fromPath(file(params.star_genome_lib, checkIfExists: true))
            .collect()
            .set { star_genome_lib }
    }
    // Prepare chimeric junctions files and input into STAR-fusion
    STAR_PREP_FUSION(star_genome_lib, fqs_ch)
    // Run STAR-fusion detection
    STAR_FUSION(star_genome_lib, fqs_ch, STAR_PREP_FUSION.out.chimera)
    md5_star(STAR_PREP_FUSION.out.bam)

    // CICERO requires GRCh37-lite or GRCh38_no_alt aligned BAMs,
    // STAR-aligner must be re-run on these specific assemblies 
    STAR_ALIGNER(fqs_ch, star_index)
    SAMTOOLS_INDEX(STAR_ALIGNER.out.bam)
    md5_cicero(STAR_ALIGNER.out.bam)
    STAR_ALIGNER.out.bam
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
    FASTQC.out.fastqc
        .concat(STAR_ALIGNER.out.log)
        .concat(STAR_PREP_FUSION.out.log)
        .collect()
        .set { mqc_ch }
    mqc_ch.view {"the mqc channel is $it"}
    MULTIQC(mqc_ch, sample_sheet)
}