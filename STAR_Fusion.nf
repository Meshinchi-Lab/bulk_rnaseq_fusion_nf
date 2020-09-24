#!/usr/bin/env nextflow

// Define the input paired fastq files in a sample sheet and genome references
// sample_sheet is tab separated with column names "Sample","R1","R2"
fqs_ch = Channel.fromPath(file(params.sample_sheet))
						.splitCsv(header: true, sep: '\t')
						.map { sample -> [sample["Sample"] + "_", file(sample["R1"]), file(sample["R2"])]}
genome_lib = params.genome_lib


// define the output directory .
params.output_folder = "./starfusion/"


//Run star-fusion on all fastq pairs and save output with the sample ID
process STAR_Fusion {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "trinityctat/starfusion:1.8.1"
	cpus 16
	memory "126 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path genome_lib
	tuple val(Sample), file(R1), file(R2) from fqs_ch

	//define output files to save to the output_folder by publishDir command
	output:
	file "*Aligned.sortedByCoord.out.bam"
	file "*SJ.out.tab"
	file "*Log.final.out"
	file "*Chimeric.out.junction"
	file "${Sample}/*abridged.coding_effect.tsv"
	path "${Sample}/FusionInspector-inspect" optional true

	"""
	set -eou pipefail

	STAR --genomeDir \$PWD/$genome_lib/ref_genome.fa.star.idx \
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


	/usr/local/src/STAR-Fusion/STAR-Fusion --genome_lib_dir \$PWD/$genome_lib \
	  	--chimeric_junction "${Sample}Chimeric.out.junction" \
			--left_fq $R1 \
			--right_fq $R2 \
	  	--CPU 8 \
	  	--FusionInspector inspect \
	  	--examine_coding_effect \
	  	--denovo_reconstruct \
	  	--output_dir $Sample

	#make dummy output directory+files, since samples that have no fusion calls will not make an output directory from Fusion Inspector
	#mkdir -p ${Sample}/FusionInspector-inspect
	#touch ${Sample}/FusionInspector-inspect/file{1..3}.txt

	echo ----------------------------------------
	echo "list all output files in $PWD"
	ls -1 \$PWD

	echo -----------------------------------------
	echo "list all output files in sample directory"
	ls -1 $Sample

	"""
}
