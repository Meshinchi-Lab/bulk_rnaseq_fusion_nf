#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s-eco-public"


# Load the module
# ml nextflow/19.11.0-edge
ml nextflow/20.04.0-edge

#I need to change some parameters so that the BAM file is an optional output with a flag, since I will not often really need it.

#More Info:
#https://github.com/FusionInspector/FusionInspector/wiki/FusionInspector-Outputs-Described

#Execute the next flow workflow  $BASE_BUCKET/work
nextflow run -c ~/nextflow.config STAR_Fusion.nf \
    --sample_sheet sample_sheets/ \
    --genome_lib $BASE_BUCKET/Reference_Data/XXXXXX/ctat_genome_lib_build_dir/ \
    --output_folder  $BASE_BUCKET/CSU_Canine_AML/RNAseq_Illumina_Data/StarFusion \
    -with-report STARFusion_CSU_Canine_AML_report.html \
    -cache  TRUE \
    -resume
