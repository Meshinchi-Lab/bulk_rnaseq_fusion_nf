#!/usr/bin/env nextflow

// Define the input genome reference files in a csv file
// inputs.txt is tab separated with column names "GENOME","GTF","DFAM","PFAM"
//the file must include full file paths (bucket and prefix) of thier location on S3.
inputs_ch = Channel.fromPath(file(params.inputs))
						.splitCsv(header: true, sep: '\t')
						.map { refs -> [file(refs["GENOME"]), file(refs["GTF"]), file(refs["DFAM"]), file(refs["PFAM"])]}


// define the output directory .
params.output_folder = "./STAR_genome_builds/"


//build a CTAT resource library for STAR-Fusion use.
process build_genome_refs {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "trinityctat/starfusion:1.10.0"
	cpus 16
	memory "126 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	tuple file(GENOME), file(GTF), file(DFAM), file(PFAM) from inputs_ch

	//define output files to save to the output_folder by publishDir command
	output:
	path "ctat_genome_lib_build_dir"

	"""
	set -eou pipefail
	echo \$STAR_FUSION_HOME

	\$STAR_FUSION_HOME/ctat-genome-lib-builder/prep_genome_lib.pl \
	                       --genome_fa $GENOME \
	                       --gtf $GTF \
												 --dfam_db $DFAM \
	                       --pfam_db $PFAM \
	                       --output_dir \$PWD
	                       --CPU 16

	"""
}
