#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s/SR"


# Load the module
# ml nextflow/19.11.0-edge
ml nextflow/20.04.0-edge


#Execute the next flow workflow  $BASE_BUCKET/work
nextflow run -c ~/nextflow.config CICERO_Fusion.nf \
    --sample_sheet sample_sheets/ \
    --genome_lib $BASE_BUCKET/Reference_Data/CICERO_Ref/  \
    --output_folder  $BASE_BUCKET/starfusion/SJ_Data/ \
    -with-report STAR-Fusion_SJ_RNAseq_report.html \
    -work-dir $BASE_BUCKET/work \
    -cache  TRUE \
    -resume
