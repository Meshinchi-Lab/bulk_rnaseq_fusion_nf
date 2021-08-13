#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s-eco-public"


# Load the module
# ml nextflow/19.11.0-edge
ml nextflow/20.04.0-edge

#docker run  -ti --rm --name rseqc  --entrypoint /bin/bash quay.io/biocontainers/rseqc:4.0.0--py39h38f01e4_1
#Execute the nextflow workflow
nextflow run -c ~/nextflow.config rseqc_TIN.nf \
    --sample_sheet sample_sheets/CellLines_Sample_Sheet.txt \
    --gene_model $BASE_BUCKET/Reference_Data/RSeQC_Ref/hg19_Ensembl_gene.bed \
    --output_folder  $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/RSeQC/ \
    -with-report mRNA_TIN_report.html \
    -cache  TRUE \
    -resume
