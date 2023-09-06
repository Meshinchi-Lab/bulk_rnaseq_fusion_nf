// Run fastQC to check each input fastq for quality metrics
// from nextflow training April 2020 at Fred Hutch, Seattle WA
process fastqc {

    //use image on quay.io
    container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"

    // if process fails, retry running it
    errorStrategy "retry"

    input:
    tuple val(sample), file(R1), file(R2)

    output:
    path "${sample}", emit: fastqc

    script:
    """
    #mkdir fastqc_${sample}
    fastqc -o ${sample} -t ${task.cpus} -f fastq -q $R1 $R2
    """
}


//Run multiQC to concatenate the results of the fastQC process
process multiqc {

    //use image on quay.io
    // container "quay.io/lifebitai/multiqc:latest"
    container "quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0"

    // if process fails, retry running it
    errorStrategy "retry"

    input:
    path "*"
    val sample_sheet

    output:
    path "${sample_sheet}_multiqc_report.html"

    script:
    """
    set -eou pipefail

    multiqc -v --filename "${sample_sheet}_multiqc_report.html"  .
    """
}


//Run star-fusion on all fastq pairs and save output with the sample ID
process STAR_Prep_Fusion {

    // use TrinityCTAT image repo on Quay.io from Biocontainers
    container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"

    // declare the input types and its variable names
    input:
    path genome_lib
    tuple val(sample), file(R1), file(R2)

    //define output files to save to the output_folder by publishDir command
    output:
    tuple val(sample), path("*Aligned.sortedByCoord.out.bam")    , emit: bam
    path("*SJ.out.tab")                                          , emit: juncs
    path("*Log.final.out")                                       , emit: log
    path("*Chimeric.out.junction")                               , emit: chimera

    script:
    """
    set -eou pipefail

    STAR --runMode alignReads \
        --genomeDir "${genome_lib}/ref_genome.fa.star.idx" \
        --runThreadN ${task.cpus} \
        --readFilesIn $R1 $R2 \
        --outFileNamePrefix "${sample}" \
        --outReadsUnmapped None \
        --twopassMode Basic \
        --twopass1readsN -1 \
        --readFilesCommand "gunzip -c" \
        --outSAMunmapped Within \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI NM MD AS nM jM jI XS \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
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
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --chimFilter banGenomicN
    """
}

//Run star-fusion on all fastq pairs and save output with the sample ID
process STAR_Fusion {
    // use TrinityCTAT image repo on Quay.io from Biocontainers
    container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"

    // declare the input types and its variable names
    input:
    path genome_lib
    tuple val(sample), file(R1), file(R2)
    path chimeric_juncs

    //define output files to save to the output_folder by publishDir command
    output:
    path("${sample}/*abridged.coding_effect.tsv")   , emit: fusions
    path("${sample}/FusionInspector-inspect")       , emit: inspector, optional: true

    script:
    """
    set -eou pipefail
    #--denovo_reconstruct
    STAR-Fusion --genome_lib_dir \$PWD/$genome_lib \
        --chimeric_junction "${chimeric_juncs}" \
        --left_fq $R1 \
        --right_fq $R2 \
        --CPU ${task.cpus} \
        --tmpdir "\$PWD" \
        --extract_fusion_reads \
        --FusionInspector inspect \
        --examine_coding_effect \
        --output_dir $sample \
        --verbose_level 2
    """
}

//build a CTAT resource library for STAR-Fusion use.
process build_genome_refs {

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

//Build GRCh37-lite index for CICERO 
process STAR_index {

    // use image on quay.io
    container "quay.io/jennylsmith/starfusion:1.8.1"

    // if process fails, retry running it
    errorStrategy "retry"

    //input genome fasta and gtf
    input: 
    path fasta
    path gtf
    
    //output the index into a diretory, and the logfile
    output:
    path "GenomeDir", emit: index
    path "Log.out", emit: log

    script:
    """
    set -eou 

    mkdir \$PWD/GenomeDir
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir \$PWD/GenomeDir \
        --genomeFastaFiles $fasta \
        --sjdbGTFfile $gtf
    """
}

process STAR_aligner {

    // use TrinityCTAT image repo on Quay.io from Biocontainers
    container "quay.io/jennylsmith/starfusion:1.8.1"
    label 'star_increasing_mem'

    input:
    path star_index_out
    tuple val(sample), file(R1), file(R2)

    output:
    tuple val(sample), path("*.bam")            , emit: bam
    path "*SJ.out.tab"                          , emit: juncs
    path "*Log.final.out"                       , emit: log

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // def seq_center      = seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$seq_center' 'SM:$prefix' $seq_platform " : "--outSAMattrRGline ID:$prefix 'SM:$prefix' $seq_platform "
    """
    set -eou pipefail 

    STAR --runMode alignReads \
        --genomeDir  "\$PWD/${star_index_out}/GenomeDir" \
        --runThreadN ${task.cpus} \
        --readFilesIn $R1 $R2 \
        --outFileNamePrefix ${sample} \
        --outReadsUnmapped None \
        --twopassMode Basic \
        --twopass1readsN -1 \
        --readFilesCommand "gunzip -c" \
        --outSAMunmapped Within \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI NM MD AS nM jM jI XS 
    """
}

//Run CICERO fusion detection on all bam files and save output with the sample ID
process CICERO {

    // use CICERO repo on docker hub.
    // container "quay.io/jennylsmith/cicero:df59166"
    container "ghcr.io/stjude/cicero:v1.9.6"

    // if process fails, retry running it
    errorStrategy "retry"

    // declare the input types and its variable names
    input:
    tuple val(sample), file(BAM), file(BAI)
    path cicero_genome_lib
    val genome

    //define output files to save to the output_folder by publishDir command
    output:
    // path "${sample}/CICERO_DATADIR/*/*.txt" optional true
    path("${sample}/*/*final_fusions.txt"), emit: cicero
    path("${sample}/*/*.txt")
    path("${sample}/*.log")

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -eou pipefail
    export TMPDIR="\$PWD"

    # CICERO fusion detection algorithm
    Cicero.sh -n ${task.cpus} \
        -b \$PWD/$BAM \
        -g ${genome} \
        -r \$PWD/$cicero_genome_lib \
        -o ${sample}
    """
}

//Helper function to unzip files when needed 
process unzip {

  input:
    path zipped_file  

    output:
    path "*", emit: unzipped_file

    script:
    """
    gunzip -f $zipped_file 
    """
}

//Create MD5sum checks for all files in the channel
process MD5sums {

    container "centos/centos:centos7"

    // if process fails, retry running it
    errorStrategy "retry"

    // declare the input types and its variable names
    input:
    tuple val(meta), path(input)

    //define output files to save to the output_folder by publishDir command
    output:
    path "*.md5"

    script:
    """
    set -eou pipefail
    echo "Creating MD5sum checks"
    md5sum ${input} > ${input}.md5
    """
}
