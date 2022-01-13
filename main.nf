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


include { STAR_Fusion; fastqc; multiqc; CICERO } from './fusion-processes.nf'

workflow  {

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

    //CICERO needs to be dependent on the STAR_Fusion.out so its done sequentially.
    STAR_Fusion.out.BAM
        .map { BAM  -> [BAM.baseName.split(/_Aligned.+/)[0], file(BAM)] }
        .set { bam_ch }

    //Run CICERO on the STAR-aligner BAM files.
    //CICERO(params.cicero_genome_lib, bam_ch)
}
