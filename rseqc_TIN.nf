#!/usr/bin/env nextflow

// Define the input BAM files, sorted and indexed, in a sample sheet and the genome references
// sample_sheet is tab separated with column names "Sample","BAM", "INDEX"
bam_ch = Channel.fromPath(file(params.sample_sheet))
						.splitCsv(header: true, sep: '\t')
						.map { inputs -> [inputs["Sample"], file(inputs["BAM"]), file(inputs["INDEX"])] }

//Gene Model references
gene_model = params.gene_model

// define the output directory .
params.output_folder = "./rseqc/"


//Run tin.py on all BAM files and save output with the sample ID
process tin_scores {

	publishDir "$params.output_folder/"

	// use Bioconainers repo on quay.io.
	container "quay.io/biocontainers/rseqc:4.0.0--py39h38f01e4_1"
	cpus 4
	memory "16 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path gene_model
	tuple val(Sample), file(BAM), file(INDEX) from bam_ch

	//define output files to save to the output_folder by publishDir command
	output:
	file "*"

	"""
	set -eou pipefail
	ls -1 \$PWD

	tin.py -i $BAM -r $gene_model
	rm $BAM #to avoid upload to workDir

	echo ----------------------------------------
	ls -1 \$PWD

	"""
}
