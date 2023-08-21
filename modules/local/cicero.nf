process CICERO {

	publishDir "$params.CICERO_out"

	// use CICERO repo on docker hub.
	container "quay.io/jennylsmith/cicero:df59166"
	cpus 2
	memory "16 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path cicero_genome_lib
	tuple val(Sample), file(BAM)

	//define output files to save to the output_folder by publishDir command
	output:
	path "${Sample}/CICERO_DATADIR/*/*.txt" optional true

 	script:
	"""
	set -eou pipefail

	#list all files in the container
	echo  -------------
	echo "the bam file is" $BAM
	ls -alh
	echo  -------------

	#index the bam file
	samtools index $BAM

	#run CICERO fusion detection algorithm
	Cicero.sh -n 2 -b $BAM -g "GRCh37-lite" \
			-r \$PWD/$cicero_genome_lib/ \
			-o ${Sample}

	echo -----------------------------------------
	echo "list all output files in sample directory"
	ls -1d $Sample/CICERO_DATADIR/
	"""
}
