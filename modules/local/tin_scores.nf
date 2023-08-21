process TIN_SCORES {

	publishDir "$params.TIN_SCORES/"

	// use Bioconainers repo on quay.io.
	container "quay.io/biocontainers/rseqc:4.0.0--py39h38f01e4_1"
	cpus 2
	memory "16 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
  //Sample == SAMPLE_ID, BAM == location of BAM file, BAI == location of BAM file .bai index
  //gene_model is pre-made reference file found at gene_model=$BASE_BUCKET/Reference_Data/RSeQC_Ref/hg19_Ensembl_gene.bed
	input:
	path gene_model
	tuple val(Sample), file(BAM), file(BAI)

	//define output files to save to the output_folder by publishDir command
	output:
	file "*"

	script:
	"""
	set -eou pipefail

	#Compute transcript integrity scores
		tin.py -i $BAM -r $gene_model
	#to avoid uploading BAM to workDir
	rm $BAM
	"""
}
