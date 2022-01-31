nextflow.enable.dsl=2


// define the output directories.
params.STAR_Fusion = "./starfusion/"
params.CICERO_out = "./CICERO/"
params.multiQC_out = "./multiqc/"


//Define message for the process logs.
log.info """\
         R N A S E Q - F U S I O N  P I P E L I N E
         ===================================
         transcriptome: ${params.star_genome_lib}
         transcriptome: ${params.cicero_genome_lib}
         samples      : ${params.sample_sheet}
         """
         .stripIndent()


include { STAR_Fusion; fastqc; multiqc; STAR_index; STAR_aligner; CICERO } from './fusion-processes.nf'


workflow star_index {
    Channel.fromPath(params.gtf_url).set{gtf}
    STAR_index(path(params.genome_fasta), gtf)
}

//take: is defining the name space 
workflow star_fusion {

    main:
    // Define the input paired fastq files in a sample sheet and genome references.
    //The sample_sheet is tab separated with column names "Sample","R1","R2"
    fqs_ch = Channel.fromPath(file(params.sample_sheet))
                    .splitCsv(header: true, sep: '\t')
                    .map { sample -> [sample["Sample"] + "_", file(sample["R1"]), file(sample["R2"])]}
    //processes are treated like functions
    STAR_Fusion(params.star_genome_lib, fqs_ch)
    //run QC on the fastq files
    fastqc(fqs_ch)
    sample_sheet=file(params.sample_sheet)
    multiqc(fastqc.out.collect(), sample_sheet.simpleName)
}

workflow  cicero {
    main: 
    // Define the input paired fastq files in a sample sheet and genome references.
    //The sample_sheet is tab separated with column names "Sample","R1","R2"
    fqs_ch = Channel.fromPath(file(params.sample_sheet))
                    .splitCsv(header: true, sep: '\t')
                    .map { sample -> [sample["Sample"] + "_", file(sample["R1"]), file(sample["R2"])]}
    // Align the fastq files to GRCh37-lite
    STAR_aligner(params.star_index_out, fqs_ch)
    //CICERO requires GRCh37-lite aligned BAMs, so dependent on STAR-aligner BAM 
    STAR_aligner.out.BAM
        .map { BAM  -> [BAM.baseName.split(/_Aligned.+/)[0], file(BAM)] }
        .set { bam_ch }
    //Run CICERO on the STAR-aligner BAM files.
    CICERO(params.cicero_genome_lib, bam_ch)
}