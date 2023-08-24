#!/bin/bash

set -eu
DATE=$(date +%F)
NFX_CONFIG=./nextflow.config
#Options: 'local_apptainer', 'PBS_apptainer', 'local_singularity', 'PBS_singularity'
NFX_PROFILE='PBS_apptainer'
#Options: star_index, fusion_calls
NFX_ENTRY='fusion_calls'
#The output prefix on filenames for reports/logs
REPORT=${1:-"pipeline_report"}

# Load the modules 
if [[ $NFX_PROFILE =~ "singularity" ]]
then
    module load singularity
elif [[ $NFX_PROFILE =~ "apptainer" ]]
then
    module load apptainer
fi

#More Info:
#https://github.com/FusionInspector/FusionInspector/wiki/FusionInspector-Outputs-Described

# Execute the nextflow index workflow
PREFIX=${REPORT}_${DATE}
nextflow -c ${NFX_CONFIG}\
    -log reports/${PREFIX}_nextflow.log \
    run main.nf \
    -entry ${NFX_ENTRY} \
    -profile ${NFX_PROFILE} \
    -with-report reports/${PREFIX}.html \
    -with-dag dag/${PREFIX}_dag.pdf \
    -cache TRUE \
    -resume
