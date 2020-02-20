#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s/SR"


# Load the module
ml nextflow/19.11.0-edge


#Execute the next flow workflow  $BASE_BUCKET/work
NXF_VER=19.10.0
nextflow run -c ~/nextflow.config STAR_Fusion.nf \
    --sample_sheet sample_sheets/relapse_sample_sheet.txt \
    --genome_lib $BASE_BUCKET/starfusion/GRCh37_gencode_v19_CTAT_lib_Oct012019/ctat_genome_lib_build_dir \
    --output_folder  $BASE_BUCKET/starfusion/relapsed_AML \
    -with-report STAR-Fusion_relapse_report.html \
    -work-dir $BASE_BUCKET/work \
    -cache  TRUE \
    -resume
