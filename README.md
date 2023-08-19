# STAR-fusion NF

 STAR-fusion run using the Nextflow workflow manager.

 This workflow requires:
  1. Sample sheet, tab delimited
  2. The location of the CTAT resource library
  3. The location of CICERO genome resource library

The output of STAR aligner, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki), and [Fusion-Inspector](https://github.com/FusionInspector/FusionInspector/wiki) will be uploaded to an S3 bucket. This includes the most relevant output files, such as SJ.out.tab, aligned.bam, and chimeric.junctions.tab, and the fusion inspector HTML report. In addition, the fastq files undergo quality control checks and a multiQC report is generated and uploaded an S3 bucket.  The workflow also includes the [CICERO](https://github.com/stjude/CICERO) fusion detection algorithm that is run using the aligned.bam from STAR-aligner output.  

The output files will be put into a directory that is named after the sample ID provided in the sample sheet file.  

The STAR-Fusion docker image can be updated easily by selecting the latest image from either 1) the CTAT [docker hub](https://hub.docker.com/r/trinityctat/starfusion) or [quay.io](quay.io). Then update the `fusion-processes.nf` with the appropriate image and tag. The containers utilized here were developed by `biocontainers` repository on quay.io, and images *must not* include `ENTRYPOINT` or it may cause errors when executed through AWS S3.  

The pre-build STAR-Fusion genome references can be found [here](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/__genome_libs_StarFv1.10/)

# To Run

## First Step: Sample Sheet

First, create a sample manifest for the fastq files that are hosted in an S3 bucket. The manifest file is a simple 3 column, tab-delimited text file with the column names "Sample", "R1", and "R2", where R1 and R2 are paired end-fastqs.

I have created a sample manifest file that can be used to select the appropriate files, with a demonstration in `Nextflow_AWS_Sample_Sheets_from_Manifest.ipynb` that uses the associated `create_sample_sheet.py`. Below is a very simple example directly on the command line.

```
python3 create_sample_sheet.py "my-s3-bucket" "Fastq" --prefix "TARGET_AML/RNAseq_Illumina_Data/Fastq/"
```

This can also be accomplished using R with 

```
library(aws.s3)
library(aws.signature)
library(tidyr) 

creds <- aws.signature::use_credentials(profile = "default")
Sys.setenv("AWS_ACCESS_KEY_ID" = creds$default$AWS_ACCESS_KEY_ID,
           "AWS_SECRET_ACCESS_KEY" = creds$default$AWS_SECRET_ACCESS_KEY,
           "AWS_DEFAULT_REGION"="us-west-2")
BUCKET="my-s3-bucket"
PREFIX="TARGET_AML/RNAseq_Illumina_Data/Fastq"
fastqs <- get_bucket_df(bucket = BUCKET, 
                        prefix = PREFIX,
                        max = Inf) %>%
           mutate(filename=str_split_fixed(Key, pattern = "/", n=4)[,4],
                  Sample=basename(filename),
                  Read=case_when(
                       grepl("_r[1].fq.gz|_R[1]_.+|r[1].fastq.gz", filename) ~ "R1", 
                       grepl("_r[2].fq.gz|_R[2]_.+|r[2].fastq.gz", filename) ~ "R2")) %>% 
           pivot_wider(id_cols=c(Sample), 
                 names_from=Read, 
                 values_from=filename)                      
```


## Second Step: Execute the Workflow

Edit the `main_run.sh` file to contain the correct output directories and point to the sample sheet in `--sample_sheet`. You can also update the filename for the output html report in `-with-report`. Then  simply call the script on the command line.

The nextflow config file can be specified by upating the variable `NFX_CONFIG` for the workflow to be run locally, on AWS Batch, or using Singularity on an HPC system. 

note: I have access to an HPC with specific software modules that can be loaded. If you have a custom installation of Nextflow, make sure the nextflow executable is in your PATH.

```
./main_run.sh
```

# CICERO Fusion Detection

This repo also provides a NF workflow for the [CICERO](https://github.com/stjude/CICERO) v1.7.1 fusion detection algorithm. The image was build from the Dockerfile provided on the CICERO github with minimal modification by removing the `ENTRYPOINT` command in order to be compatible with AWS Batch.

The genome references for CICERO can be found [here](https://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa) and [here](https://doi.org/10.5281/zenodo.3817656).

Docker Image: https://hub.docker.com/repository/docker/jennylsmith/cicero 

## Quality Control 

In addition, running [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiQC](https://multiqc.info/) on the input fastq files. 

# DSL2 Resources

* [Tutorial](https://github.com/nextflow-io/nfcamp-tutorial)
* [Documentation](https://www.nextflow.io/docs/latest/dsl2.html)
* [Dynamic memory allocation](https://lucacozzuto.medium.com/handling-failing-jobs-with-nextflow-24405b97b679)

# Resources on Automated Builds
* [Fred Hutch SciWiki](https://sciwiki.fredhutch.org/hdc/hdc_building_containers/#step-2-define-your-docker-image)