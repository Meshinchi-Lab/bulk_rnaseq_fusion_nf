#!/usr/bin/env nextflow

// Define the input paired fastq files in a sample sheet and genome references
// sample_sheet is tab separated with column names "Sample","Filename"
input_ch = Channel.fromPath(file(params.sample_sheet))
						.splitCsv(header: true, sep: '\t')
						.map { inputs -> [inputs["Sample"], file(inputs["Filename"]) ]}


// define the output directory .
params.output_folder = "./MD5sums/"


//Run star-fusion on all fastq pairs and save output with the sample ID
process MD5sums {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "ubuntu:latest"
	cpus 4
	memory "16 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	tuple val(Sample), file(Filename) from input_ch

	//define output files to save to the output_folder by publishDir command
	output:
	file "*.md5"

	"""
	set -eou pipefail
	ls -1 \$PWD

	echo "Creating MD5sum checks"
	hashes=\$(echo $Filename)
	hashes=\${hashes}.md5
	md5sum $Filename > \$hashes


	echo ----------------------------------------
	ls -1 \$PWD

	"""
}
