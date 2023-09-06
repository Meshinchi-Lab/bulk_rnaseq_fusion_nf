process BUILD_GENOME_REFS {

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
