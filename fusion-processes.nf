// Run fastQC to check each input fastq for quality metrics
// from nextflow training April 2020 at Fred Hutch, Seattle WA
process fastqc {

    //use image on quay.io
    container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
    cpus 4
    memory "32 GB"

    // if process fails, retry running it
    errorStrategy "retry"

    input:
    tuple val(Sample), file(R1), file(R2)

    output:
    path "fastqc_${Sample}_logs"

    script:
    """
    mkdir fastqc_${Sample}_logs
    fastqc -o fastqc_${Sample}_logs -t 6 -f fastq -q $R1 $R2
    """
}


//Run multiQC to concatenate the results of the fastQC process
// from nextflow training April 2020 at Fred Hutch, Seattle WA
// mode:'copy'
process multiqc {

    publishDir params.multiQC

    //use image on quay.io
    container "quay.io/lifebitai/multiqc:latest"
    cpus 4
    memory "32 GB"

    // if process fails, retry running it
    errorStrategy "retry"

    input:
    path "*"

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc -v .
    """
}


//Run star-fusion on all fastq pairs and save output with the sample ID
process STAR_Fusion {

	publishDir "$params.output_folder/"

	// use TrinityCTAT repo on docker hub.
	container "trinityctat/starfusion:1.8.1"
	cpus 8
	memory "64 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path genome_lib
	tuple val(Sample), file(R1), file(R2)

	//define output files to save to the output_folder by publishDir command
	output:
	file "*Aligned.sortedByCoord.out.bam"
	file "*SJ.out.tab"
	file "*Log.final.out"
	file "*Chimeric.out.junction"
	file "${Sample}/*abridged.coding_effect.tsv"
	path "${Sample}/FusionInspector-inspect" optional true

  script:
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


	echo ----------------------------------------
	echo "list all output files in $PWD"
	ls -1 \$PWD

	echo -----------------------------------------
	echo "list all output files in sample directory"
	ls -1 $Sample

	"""
}


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
	tuple file(GENOME), file(GTF), val(DFAM), val(PFAM)

	//define output files to save to the output_folder by publishDir command
	output:
	path "ctat_genome_lib_build_dir"

  script:
	"""
	set -eou pipefail
	echo \$STAR_FUSION_HOME
	ls -alh \$PWD

	\$STAR_FUSION_HOME/ctat-genome-lib-builder/prep_genome_lib.pl \
	                       	--genome_fa \$PWD/$GENOME \
			 	                  --gtf \$PWD/$GTF \
				                  --dfam_db $DFAM \
	                       	--pfam_db $PFAM \
	                       	--CPU 16

	find . -name "ctat_genome_lib_build_dir" -type d

	"""
}


//Run CICERO fusion detection on all bam files and save output with the sample ID
process CICERO {

	publishDir "$params.CICERO"

	// use CICERO repo on docker hub.
	container "jennylsmith/cicero:v1.7.1"
	cpus 8
	memory "64 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path cicero_genome_lib
	tuple val(Sample), file(BAM)

	//define output files to save to the output_folder by publishDir command
	output:
	path "${Sample}" optional true

  script:
	"""
	set -eou pipefail

  #index the bam file
  samtools index $BAM

  #run CICERO fusion detection algorithm
	Cicero.sh -n 8 -b \$PWD/$BAM -g "GRCh37-lite" \
		-r \$PWD/$cicero_genome_lib/ \
		-o ${Sample}

	echo ----------------------------------------
	echo "list all output files in $PWD"
	ls -1 \$PWD

	echo -----------------------------------------
	echo "list all output files in sample directory"
	ls -1 $Sample
	"""
}

//Run tin.py for QC check on all BAM files and save output with the sample ID
process tin_scores {

	publishDir "$params.R/"

	// use Bioconainers repo on quay.io.
	container "quay.io/biocontainers/rseqc:4.0.0--py39h38f01e4_1"
	cpus 4
	memory "16 GB"

	// if process fails, retry running it
	errorStrategy "retry"

	// declare the input types and its variable names
	input:
	path gene_model
	tuple val(Sample), file(BAM), file(INDEX)

	//define output files to save to the output_folder by publishDir command
	output:
	file "*"

  script:
	"""
	set -eou pipefail
	ls -1 \$PWD

	tin.py -i $BAM -r $gene_model
	rm $BAM #to avoid upload to workDir

	echo ----------------------------------------
	ls -1 \$PWD
	"""
}



//Create MD5sum checks for all files in the channel
process MD5sums {

	// use ubuntu repo on docker hub.
	container "ubuntu:latest"
	cpus 4
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
