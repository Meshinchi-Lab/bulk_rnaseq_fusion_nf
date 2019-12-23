#!/usr/bin/env nextflow

// Define the input paired fastq files in a sample sheet and genome references
// sample_sheet is tab separated with column names "Sample","R1","R2"
fqs_ch = Channel.fromPath(file(params.sample_sheet))
						.splitCsv(header: true, sep: '\t')
						.map { sample -> [sample["Sample"], file(sample["R1"]), file(sample["R2"])]}
genome_lib = params.genome_lib

// define the output directory .
params.output_folder = "./starfusion/"

//Run star-fusion on all fastq pairs and save output with the sample ID
process STAR_Fusion {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "trinityctat/starfusion:latest"
	cpus 8
	memory "64 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path genome_lib
	tuple val(Sample), file(R1), file(R2) from fqs_ch

	//define output files to save to the output_folder by publishDir command
	output:
	file "*Log.final.out"
	file "*Log.out"
	file "*Aligned.out.bam"
	file "*SJ.out.tab"
	file "*Chimeric.out.junction"
	file "*ReadsPerGene.out.tab"
	file "*.abridged.coding_effect.tsv"
	path "*FusionInspector*" maxDepth 1

	"""
	set -eou pipefail
	echo \$PWD/$genome_lib

	/usr/local/src/STAR-2.7.2b/bin/Linux_x86_64/STAR --runMode alignReads \
		--runThreadN 8 \
		--genomeDir \$PWD/$genome_lib  \
		--readFilesIn $R1 $R2 \
		--readFilesCommand zcat \
		--outReadsUnmapped None \
		--outSAMunmapped Within  \
		--twopassMode Basic \
		--twopass1readsN -1 \
		--outFileNamePrefix "${R1.simpleName}" \
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


	/usr/local/src/STAR-Fusion/STAR-Fusion \
		--genome_lib_dir \$PWD/$genome_lib  \
		-J Chimeric.out.junction
		--CPU 8 \
		--FusionInspector inspect \
		--examine_coding_effect \
		--denovo_reconstruct \
		--output_dir "${R1.simpleName}"

	"""
}
