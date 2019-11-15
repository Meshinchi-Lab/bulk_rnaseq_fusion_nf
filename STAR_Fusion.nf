#!/usr/bin/env nextflow

// Define the input paired fastq files in a sample sheet and genome references
// sample_sheet is tab separated with column names "Sample","R1","R2"
fqs_ch = Channel.fromPath(file(params.sample_sheet))
						.splitCsv(header: true, sep: '\t')
						.map { sample -> [sample["Sample"], file(sample["R1"]), file(sample["R2"])]}
genome_lib = params.genome_lib_dir

// define the output directory .
params.output_folder = "./starfusion/"

//Run star-fusion on all fastq pairs and save output with the sample ID
process STAR_Fusion {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "trinityctat/starfusion:1.8.1"
	cpus 8
	memory "60 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path genome_lib
	tuple val(Sample), file(R1), file(R2) from fqs_ch

	//define output files to save to the output_folder by publishDir command
	output:
	file "*"


	"""
	set -eou pipefail

	/usr/local/src/STAR-Fusion/STAR-Fusion --left_fq $R1 --right_fq $R2 \
		--genome_lib_dir $genome_lib  \
		--CPU 8 \
		--FusionInspector validate \
		--examine_coding_effect
	"""
}
