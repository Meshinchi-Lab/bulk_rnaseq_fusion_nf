#!/bin/bash

# Nextflow Version
NXF_VER=21.04.1
NFX_CONFIG=nextflow.local.config

nextflow run -c ${NFX_CONFIG} gtf_url_test.nf