process MD5SUMS {

	// use ubuntu repo on docker hub.
	container "ubuntu:latest"
	cpus 2
	memory "16 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	file Filename

	//define output files to save to the output_folder by publishDir command
	output:
	file "*.md5"

	script:
	"""
	set -eou pipefail

	echo "Creating MD5sum checks"
	hashes=${Filename}.md5
	md5sum $Filename > \$hashes
	"""
}
