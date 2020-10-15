#!/usr/bin/env nextflow

// Define the input of BAM files in a sample sheet and genome references
// sample_sheet is tab separated with column names "Sample" and "BAM"
bam_ch = Channel.fromPath(file(params.sample_sheet))
						.splitCsv(header: true, sep: '\t')
						.map { sample -> [sample["Sample"] + "_", file(sample["BAM"])]}
genome_lib = params.genome_lib


// define the output directory .
params.output_folder = "./CICERO/"


//Run CICERO fusion detection on all bam files and save output with the sample ID
process CICERO {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "jennylsmith/cicero:v0.3.0"
	cpus 8
	memory "64 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path genome_lib
	tuple val(Sample), file(BAM) from bam_ch

	//define output files to save to the output_folder by publishDir command
	//path "${Sample}/" optional true
	output:
	file "*"

	"""
	set -eou pipefail

	Cicero.sh -n 8 -b $BAM -g "GRCh37-lite" \
		-r \$PWD/$genome_lib/ \
		-o $output_folder/$Sample

	echo ----------------------------------------
	echo "list all output files in $PWD"
	ls -1 \$PWD

	echo -----------------------------------------
	echo "list all output files in sample directory"
	ls -1 $Sample

	"""
}
