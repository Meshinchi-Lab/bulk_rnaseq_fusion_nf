process BUILD_GENOME_REFS {

	publishDir "$params.Reference_Data/"

	// use TrinityCTAT image from biocontainers
	container "quay.io/biocontainers/star-fusion:1.9.1--0"
	cpus 16
	memory "126 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	tuple file(GENOME), file(GTF), val(DFAM), val(PFAM)

	//define output files to save to the output_folder by publishDir command
	output:
	path "ctat_genome_lib_build_dir"

  script:
	"""
	set -eou pipefail
	prep_genome_lib.pl \
	        --genome_fa \$PWD/$GENOME \
			--gtf \$PWD/$GTF \
			--dfam_db $DFAM \
	        --pfam_db $PFAM \
	        --CPU 16

	find . -name "ctat_genome_lib_build_dir" -type d

	"""
}
