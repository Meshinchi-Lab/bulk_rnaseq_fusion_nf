#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s/SR"


# Load the module
# ml nextflow/19.11.0-edge
ml nextflow/20.04.0-edge


#Execute the next flow workflow  $BASE_BUCKET/work
nextflow run -c ~/nextflow.config STAR_Fusion.nf \
    --sample_sheet sample_sheets/NUP98-NAPEPLD_sample_sheet.txt \
    --genome_lib $BASE_BUCKET/starfusion/GRCh37_gencode_v19_CTAT_lib_Oct012019/ctat_genome_lib_build_dir \
    --output_folder  $BASE_BUCKET/starfusion/ \
    -with-report STAR-Fusion_NUP98-NAPEPLD_report.html \
    -work-dir $BASE_BUCKET/work \
    -cache  TRUE \
    -process.queue cpu-spot-50 \ #temporary bc spot-30 is not really working
    -resume
