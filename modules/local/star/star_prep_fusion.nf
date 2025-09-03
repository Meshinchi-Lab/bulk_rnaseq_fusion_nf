process STAR_PREP_FUSION {

    // use TrinityCTAT image repo on Quay.io from Biocontainers
    // container "quay.io/biocontainers/star-fusion:1.15.1--hdfd78af_1"

    // declare the input types and its variable names
    input:
    path genome_lib
    tuple val(sample), path(R1), path(R2)

    //define output files to save to the output_folder by publishDir command
    output:
    tuple val(sample), path("*Aligned.sortedByCoord.out.bam")    , emit: bam
    tuple val(sample), path("*SJ.out.tab")                       , emit: juncs
    tuple val(sample), path("*Chimeric.out.junction")            , emit: chimera
    path("*Log.final.out")                                       , emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sample}"
    """
    STAR --runMode alignReads \
        --genomeDir ./${genome_lib}/ref_genome.fa.star.idx \
        --runThreadN ${task.cpus} \
        --readFilesIn ${R1} ${R2} \
        --outFileNamePrefix "${prefix}" \
        $args \
        --twopassMode Basic \
        --twopass1readsN -1 \
        --outSAMunmapped Within \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \

        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 \
        --outSAMattrRGline ID:GRPundef \
        --chimMultimapScoreRange 3 \
        --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 \
        --chimNonchimScoreDropMin 10 \
        --chimFilter banGenomicN \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1
    """
}
