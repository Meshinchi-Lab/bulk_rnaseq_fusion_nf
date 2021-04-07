#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s-eco-public"


# Load the module
# ml nextflow/19.11.0-edge
ml nextflow/20.04.0-edge


#Execute the nextflow workflow
nextflow run -c ~/nextflow.config create_MD5.nf \
    --sample_sheet sample_sheets/TARGET_AML_Fastqs_MD5_Sample_Sheet_v2.txt \
    --output_folder  $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/Fastq \
    -with-report Create_MD5sums_report_v2.html \
    -cache  TRUE \
    -resume
