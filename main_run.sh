#!/bin/bash
set -e
DATE=$(date +%F)
# Change NFX_CONFIG to nextflow.aws.config if needed for AWS executor. 
NFX_CONFIG=./nextflow.singularity.config

# Load the modules
ml Singularity
export PATH=$SINGULARITYROOT/bin/:$PATH

# Declare the output directories 
FH_BASE_BUCKET="s3://fh-pi-meshinchi-s-eco-public"
STAR_FUSION_OUT=$FH_BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/STAR_Fusion/
MULTIQC_OUT=$FH_BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/multiQC/
STAR_INDEX_OUT=$FH_BASE_BUCKET/Reference_Data/CICERO_Ref/GRCh37_gencode_v19_STAR_idx/
STAR_ALIGNER_OUT=$FH_BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/CICERO/
CICERO_OUT=$FH_BASE_BUCKET/TARGET_AML/RNAseq_Illumina_Data/CICERO/

#More Info:
#https://github.com/FusionInspector/FusionInspector/wiki/FusionInspector-Outputs-Described

#Execute the nextflow index workflow
#TARGET_AML_Remission_JMML_APL_MDAnderson_sample_sheet.txt
nextflow run -c ${NFX_CONFIG} main.nf \
    --sample_sheet sample_sheets/test_sample_sheet.txt \
    --star_genome_lib $FH_BASE_BUCKET/Reference_Data/GRCh37_gencode_v19_CTAT_lib_Oct012019/ctat_genome_lib_build_dir/ \
    --STAR_Fusion_out $STAR_FUSION_OUT \
    --multiQC_out $MULTIQC_OUT \
    --fasta_file $FH_BASE_BUCKET/Reference_Data/CICERO_Ref/Homo_sapiens/GRCh37-lite/FASTA/GRCh37-lite.fa \
    --gtf_url "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz" \
    --star_index_out $STAR_INDEX_OUT \
    --STAR_aligner_out $STAR_ALIGNER_OUT \
    --cicero_genome_lib $FH_BASE_BUCKET/Reference_Data/CICERO_Ref/ \
    --CICERO_out $CICERO_OUT \
    -entry star_index \
    -with-report generate_index_${DATE}_singularity.html \
    -cache  TRUE \
    -resume
