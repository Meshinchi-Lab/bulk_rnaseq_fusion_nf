#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s-eco-public"


# Load the module
ml nextflow/20.09.0-edge

#I need to change some parameters so that the BAM file is an optional output with a flag, since I will not often really need it.

#More Info:
#https://github.com/FusionInspector/FusionInspector/wiki/FusionInspector-Outputs-Described

#Execute the next flow workflow  $BASE_BUCKET/work
nextflow run -c ~/nextflow.config main.nf \
    --sample_sheet sample_sheets/test_sample_sheet.txt \
    --star_genome_lib $BASE_BUCKET/Reference_Data/GRCh37_gencode_v19_CTAT_lib_Oct012019/ctat_genome_lib_build_dir/ \
    --output_folder  $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/STAR_Fusion \
    --multiQC $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/multiQC \
    --cicero_genome_lib $BASE_BUCKET/Reference_Data/CICERO_Ref/ \
    --CICERO $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/CICERO \
    -with-report STAR-Fusion_Remission_JMML_APL_MDAnderson_AML_report.html \
    -cache  TRUE \
    -resume
