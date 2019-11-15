#!/usr/bin/env nextflow

// Define the input paired fastq files in a sample sheet and genome references
fqs_ch = Channel.fromPath(file(params.sample_sheet))
						.splitCsv(header: true, sep: '\t')
						.map { sample -> [sample["Sample"], file(sample["R1"]), file(sample["R2"])]}
//fqs_ch.subscribe{ println "${it}\n" }
genome_lib = params.genome_lib

// define the output directory .
params.output_folder = "./starfusion/"

//Run star-fusion on all fastq pairs and save output with the sample ID
process STAR_Fusion_test {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "trinityctat/starfusion:1.8.1"
	cpus 4
	memory "30 GB"

	// if process fails, retry running it
	errorStrategy "terminate"

	// declare the input types and its variable names
	input:
	path genome_lib
	tuple val(Sample), file(R1), file(R2) from fqs_ch

	//define output files to save to the output_folder by publishDir command
	output:
	stdout "*"

	"""
	set -eou pipefail
	echo $genome_lib
	ls -alh $genome_lib
	find . -name "AnnotFilterRule.pm"
	"""
}
