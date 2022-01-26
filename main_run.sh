#!/bin/bash
set -e
BASE_BUCKET="s3://fh-pi-meshinchi-s-eco-public"

# Load the module
ml nextflow/20.09.0-edge

#I need to change some parameters so that the BAM file is an optional output with a flag, since I will not often really need it.

#More Info:
#https://github.com/FusionInspector/FusionInspector/wiki/FusionInspector-Outputs-Described

#Execute the next flow workflow
#TARGET_AML_Remission_JMML_APL_MDAnderson_sample_sheet.txt
nextflow run -c ~/nextflow.config main.nf \
    --sample_sheet sample_sheets/test_sample_sheet.txt \
    --star_genome_lib $BASE_BUCKET/Reference_Data/GRCh37_gencode_v19_CTAT_lib_Oct012019/ctat_genome_lib_build_dir/ \
    --STAR_Fusion_out  $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/STAR_Fusion/ \
    --multiQC_out $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/multiQC/ \
    --STAR_aligner_out $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/CICERO/ \
    --genome_fasta $BASE_BUCKET/Reference_Data/CICERO_Ref/Homo_sapiens/GRCh37-lite/FASTA/chromosomes/ \
    --gtf_url "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz" \
    --star_index_dir $BASE_BUCKET/Reference_Data/CICERO_Ref/GRCh37_gencode_v19_STAR_idx \
    --cicero_genome_lib $BASE_BUCKET/Reference_Data/CICERO_Ref/ \
    --CICERO_out $BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/CICERO/ \
    -with-report STAR-Fusion_Remission_JMML_APL_MDAnderson_AML_report.html \
    -cache  TRUE \
    -resume
