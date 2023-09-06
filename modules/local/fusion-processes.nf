// // Run fastQC to check each input fastq for quality metrics
// process fastqc {

//     //use image on quay.io
//     container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"

//     input:
//     tuple val(sample), file(R1), file(R2)

//     output:
//     path "${sample}", emit: fastqc

//     script:
//     def args = task.ext.args ?: ''
//     """
//     set -eou pipefail
//     mkdir ${sample}
//     fastqc \\
//         --outdir ${sample} \\
//         --threads ${task.cpus} \\
//         $args \\
//         $R1 $R2
//     """
// }


// //Run multiQC to concatenate the results of the fastQC process
// process multiqc {

//     //use image on quay.io
//     container "quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0"

//     input:
//     path(input)
//     val sample_sheet

//     output:
//     path("${sample_sheet}_multiqc_report.html")     , emit: report
//     path("*report_data")                            , emit: data, optional: true

//     script:
//     def args = task.ext.args ?: ''
//     """
//     set -eou pipefail
//     multiqc \\
//         $args \\
//         --filename "${sample_sheet}_multiqc_report.html" \\
//         \$PWD
//     """
// }


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
    def args = task.ext.args ?: ''
    """
    set -eou pipefail

    STAR --runMode alignReads \
        --genomeDir "${genome_lib}/ref_genome.fa.star.idx" \
        --runThreadN ${task.cpus} \
        --readFilesIn $R1 $R2 \
        --outFileNamePrefix "${sample}" \
        $args \
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

// //Run star-fusion on all fastq pairs and save output with the sample ID
// process STAR_Fusion {
//     // use TrinityCTAT image repo on Quay.io from Biocontainers
//     container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"

//     // declare the input types and its variable names
//     input:
//     path genome_lib
//     tuple val(sample), file(R1), file(R2)
//     path chimeric_juncs

//     //define output files to save to the output_folder by publishDir command
//     output:
//     path("${sample}/*abridged.coding_effect.tsv")   , emit: fusions
//     path("${sample}/FusionInspector-inspect")       , emit: inspector, optional: true

//     script:
//     def args = task.ext.args ?: ''
//     """
//     set -eou pipefail
//     #--denovo_reconstruct
//     STAR-Fusion \\
//         --genome_lib_dir \$PWD/$genome_lib \\
//         --chimeric_junction "${chimeric_juncs}" \\
//         --left_fq $R1 \\
//         --right_fq $R2 \\
//         --CPU ${task.cpus} \\
//         $args \\
//         --tmpdir "\$PWD"
//         --output_dir ${sample}
//     """
// }

//build a CTAT resource library for STAR-Fusion use.
process build_genome_refs {

    // use TrinityCTAT image from biocontainers
    container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"

    // declare the input types and its variable names
    input:
    tuple file(GENOME), file(GTF), val(DFAM), val(PFAM)

    //define output files to save to the output_folder by publishDir command
    output:
    path("ctat_genome_lib_build_dir"), emit: genome_dir

    script:
    def args = task.ext.args ?: ''
    """
    set -eou pipefail
    prep_genome_lib.pl \\
            --genome_fa \$PWD/$GENOME \\
            --gtf \$PWD/$GTF \\
            --dfam_db $DFAM \\
            --pfam_db $PFAM \\
            $args \\
            --CPU ${task.cpus}
    """
}

//Build GRCh37-lite or GRCh38_no_alt index for CICERO 
process STAR_index {

    // use image on quay.io
    container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"

    //input genome fasta and gtf
    input: 
    path fasta
    path gtf
    
    //output the index into a diretory, and the logfile
    output:
    path("GenomeDir"), emit: index

    script:
    def args = task.ext.args ?: ''
    """
    set -eou 

    mkdir \$PWD/GenomeDir
    STAR --runThreadN ${task.cpus} \\
        --runMode genomeGenerate \\
        $args \\
        --genomeDir \$PWD/GenomeDir \\
        --genomeFastaFiles $fasta \\
        --sjdbGTFfile $gtf
    """
}

process STAR_aligner {

    // use TrinityCTAT image repo on Quay.io from Biocontainers
    container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"
    label 'star_increasing_mem'

    input:
    tuple val(sample), file(R1), file(R2)
    path star_index_dir

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
        --genomeDir  "\$PWD/${star_index_dir}" \
        --runThreadN ${task.cpus} \
        --readFilesIn $R1 $R2 \
        --outFileNamePrefix ${sample} \
        $args \\
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

    // use CICERO container from stjude
    container "ghcr.io/stjude/cicero:v1.9.6"

    // declare the input types and its variable names
    input:
    tuple val(sample), file(BAM), file(BAI)
    path cicero_genome_lib
    val genome

    //define output files
    output:
    path("${sample}/*/*final_fusions.txt")  , emit: cicero, optional: true
    path("${sample}/*/*.txt")               , emit: outfiles
    path("${sample}/*.log")                 , emit: logs

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -eou pipefail
    export TMPDIR="\$PWD"

    # CICERO fusion detection algorithm
    Cicero.sh \\
        -n ${task.cpus} \\
        -b \$PWD/$BAM \\
        -g ${genome} \\
        -r \$PWD/$cicero_genome_lib \\
        $args \\
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
    def args = task.ext.args ?: ''
    """
    gunzip $args -f $zipped_file 
    """
}

//Create MD5sum checks for all files in the channel
process MD5sums {

    container "centos/centos:centos7"

    // declare the input types and its variable names
    input:
    tuple val(meta), path(input)

    //define output files to save to the output_folder by publishDir command
    output:
    path "*.md5"

    script:
    def args = task.ext.args ?: ''
    """
    set -eou pipefail
    echo "Creating MD5sum checks"
    md5sum $args ${input} > ${input}.md5
    """
}
