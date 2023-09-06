process STAR_ALIGNER {

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
