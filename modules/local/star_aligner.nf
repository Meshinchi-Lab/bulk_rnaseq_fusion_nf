process STAR_ALIGNER {
	publishDir "$params.STAR_ALIGNER_out/"

	// use TrinityCTAT image repo on Quay.io from Biocontainers
	container "quay.io/biocontainers/star-fusion:1.9.1--0"
	label 'star_increasing_mem'
	
	input:
	path star_index_out
	tuple val(Sample), file(R1), file(R2)

	output:
	path "*.bam", emit: BAM
	path "*SJ.out.tab"
	path "*Log.final.out"

	script:
	"""
	set -eou pipefail 

	genome_idx=\$(basename ${star_index_out})
	echo \$genome_idx

	STAR --runMode alignReads \
    	--genomeDir  \$PWD/$star_index_out \
		--runThreadN 8 \
		--readFilesIn $R1 $R2 \
		--outFileNamePrefix ${Sample} \
		--outReadsUnmapped None \
		--twopassMode Basic \
		--twopass1readsN -1 \
		--readFilesCommand "gunzip -c" \
		--outSAMunmapped Within \
		--outSAMtype BAM \
		--outSAMattributes NH HI NM MD AS nM jM jI XS 
	"""
}
