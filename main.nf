nextflow.enable.dsl=2


// define the output directories.
params.output_folder = "./starfusion/"
params.CICERO = "./CICERO/"
params.multiQC = "./multiqc/"


//Define message for the process logs.
log.info """\
         R N A S E Q - F U S I O N  P I P E L I N E
         ===================================
         transcriptome: ${params.star_genome_lib}
         transcriptome: ${params.cicero_genome_lib}
         samples      : ${params.sample_sheet}
         outdir       : ${params.output_folder}
         """
         .stripIndent()


include { STAR_Fusion; MD5sums; fastqc; multiqc; CICERO } from './fusion-processes.nf'

workflow  {
    /*
		// Define the input paired fastq files in a sample sheet and genome references.
		//The sample_sheet is tab separated with column names "Sample","R1","R2"
		fqs_ch = Channel.fromPath(file(params.sample_sheet))
				.splitCsv(header: true, sep: '\t')
				.map { sample -> [sample["Sample"] + "_", file(sample["R1"]), file(sample["R2"])]}

    //flattened channel for MD5sums to create 1 file per fastq
    files_ch = Channel.fromPath(file(params.sample_sheet))
                      .splitCsv(header: true, sep: '\t')
                      .map { Filename -> [file(Filename["R1"]), file(Filename["R2"])]}
                      .flatten()

		//processes are treated like functions
    STAR_Fusion(params.star_genome_lib, fqs_ch)
    MD5sums(files_ch)


    //run QC on the fastq files
    //direcly call a process on the output of a previous task
    fastqc(fqs_ch)
    multiqc(fastqc.out.collect())
    */


    //Create new channel for the BAM output of STAR fusion
    bam_ch = bam_ch = Channel.fromPath(params.output_folder + "*", type: 'file', checkIfExists: true)
                .filter( ~/^.+bam/ )
                .map { BAM  -> [BAM.baseName.split(/_Aligned.+/)[0], file(BAM)] }

    //Run CICERO on the STAR aligner BAM files
    CICERO(params.cicero_genome_lib, bam_ch)

}
