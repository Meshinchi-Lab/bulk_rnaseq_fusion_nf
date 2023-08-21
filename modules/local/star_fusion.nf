process STAR_FUSION {

	publishDir "$params.STAR_FUSION_out/"

	// use TrinityCTAT image repo on Quay.io from Biocontainers
	container "quay.io/biocontainers/star-fusion:1.9.1--0"
	label 'star_increasing_mem'

	// declare the input types and its variable names
	input:
	path genome_lib
	tuple val(Sample), file(R1), file(R2)

	//define output files to save to the output_folder by publishDir command
	output:
	path "*Aligned.sortedByCoord.out.bam", emit: BAM
	path "*SJ.out.tab"
	path "*Log.final.out"
	path "*Chimeric.out.junction"
	path "${Sample}/*abridged.coding_effect.tsv"
	path "${Sample}/FusionInspector-inspect" optional true

  script:
	"""
	set -eou pipefail


	#list all files in the container
	echo  -------------
	echo "the genome lib is file is" $genome_lib
	ls -alh \$PWD/$genome_lib/ref_genome.fa.star.idx
	ls -alh
	echo  -------------

	STAR --runMode alignReads \
    	--genomeDir \$PWD/$genome_lib/ref_genome.fa.star.idx \
		--runThreadN 8 \
		--readFilesIn $R1 $R2 \
		--outFileNamePrefix $Sample \
		--outReadsUnmapped None \
		--twopassMode Basic \
		--twopass1readsN -1 \
		--readFilesCommand "gunzip -c" \
		--outSAMunmapped Within \
		--outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM 63004036730 \
		--outSAMattributes NH HI NM MD AS nM jM jI XS \
		--chimSegmentMin 12 \
		--chimJunctionOverhangMin 12 \
		--chimOutJunctionFormat 1 \
		--alignSJDBoverhangMin 10 \
		--alignMatesGapMax 100000 \
		--alignIntronMax 100000 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--outSAMattrRGline ID:GRPundef \
		--chimMultimapScoreRange 3 \
		--chimScoreJunctionNonGTAG -4 \
		--chimMultimapNmax 20 \
		--chimNonchimScoreDropMin 10 \
		--peOverlapNbasesMin 12 \
		--peOverlapMMp 0.1 \
		--chimFilter banGenomicN

	STAR-Fusion --genome_lib_dir \$PWD/$genome_lib \
	  	--chimeric_junction "${Sample}Chimeric.out.junction" \
		--left_fq $R1 \
		--right_fq $R2 \
	  	--CPU 8 \
	  	--FusionInspector inspect \
	  	--examine_coding_effect \
	  	--denovo_reconstruct \
	  	--output_dir $Sample

	"""
}
