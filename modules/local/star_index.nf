process STAR_INDEX {
	publishDir "$params.star_index_out"

	// use person
	container "quay.io/jennylsmith/starfusion:1.8.1"
	cpus 16
	memory "315 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	//input genome fasta and gtf
	input: 
	path fasta
	path gtf
	
	//output the index into a diretory, and the logfile
	output:
	path "*/GenomeDir"
	path "Log.out"

	script:
	"""
	set -eou 

	mkdir \$PWD/GenomeDir
	STAR --runThreadN 16 \
		--runMode genomeGenerate \
		--genomeDir \$PWD/GenomeDir \
		--genomeFastaFiles $fasta \
		--sjdbGTFfile $gtf
	"""
}
