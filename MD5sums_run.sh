#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s-eco-public"


# Load the module
# ml nextflow/19.11.0-edge
ml nextflow/20.04.0-edge


#Execute the nextflow workflow  
nextflow run -c ~/nextflow.config create_MD5.nf \
    --sample_sheet sample_sheets/test_sample_sheet.txt \
    --output_folder  $BASE_BUCKET/ \
    -with-report Create_MD5sums_report.html \
    -cache  TRUE \
    -resume
