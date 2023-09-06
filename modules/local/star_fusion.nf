process STAR_FUSION {
    // use TrinityCTAT image repo on Quay.io from Biocontainers
    container "quay.io/biocontainers/star-fusion:1.12.0--hdfd78af_1"

    // declare the input types and its variable names
    input:
    path genome_lib
    tuple val(sample), file(R1), file(R2)
    path chimeric_juncs

    //define output files to save to the output_folder by publishDir command
    output:
    path("${sample}/*abridged.tsv")                 , emit: fusions
    path("${sample}/*coding_effect.tsv")            , emit: coding_effect, optional: true
    path("${sample}/fusion_evidence_reads_*.fq")    , emit: fastqs, optional: true
    path("${sample}/FusionInspector.log")           , emit: log, optional: true
    path("${sample}/FusionInspector-*")             , emit: inspector, optional: true

    script:
    def args = task.ext.args ?: ''
    """
    set -eou pipefail
    STAR-Fusion \\
        $args \\
        --genome_lib_dir \$PWD/$genome_lib \\
        --chimeric_junction "${chimeric_juncs}" \\
        --left_fq $R1 \\
        --right_fq $R2 \\
        --CPU ${task.cpus} \\
        --tmpdir "\$PWD" \\
        --output_dir ${sample}
    """
}
