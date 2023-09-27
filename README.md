# STAR-fusion NF

 STAR-fusion run using the Nextflow workflow manager.

 This workflow requires:
  1. Sample sheet, comma delimited
  2. The location of the CTAT resource library
  3. The location of CICERO genome resource library

The output of STAR aligner, [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki), and [Fusion-Inspector](https://github.com/FusionInspector/FusionInspector/wiki) will be uploaded to an S3 bucket. This includes the most relevant output files, such as SJ.out.tab, aligned.bam, and chimeric.junctions.tab, and the fusion inspector HTML report. In addition, the fastq files undergo quality control checks and a multiQC report is generated and uploaded an S3 bucket.  The workflow also includes the [CICERO](https://github.com/stjude/CICERO) fusion detection algorithm that is run using the aligned.bam from STAR-aligner output.  

The output files will be put into a directory that is named after the sample ID provided in the sample sheet file.  

The STAR-Fusion docker image can be updated easily by selecting the latest image from either 1) the CTAT [docker hub](https://hub.docker.com/r/trinityctat/starfusion) or [quay.io](quay.io). Then update the `fusion-processes.nf` with the appropriate image and tag.

The pre-build STAR-Fusion genome references can be found [here](https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/). 

If a new fusion resource library must be created, for example for a non-human species, there is detailed documentation available [here](https://github.com/NCIP/ctat-genome-lib-builder/wiki).

[09/26/2023]
STAR-fusion --denovo_reconstruct currently does not work. 

# To Run

## Install Nextflow

Install nextflow using the conda env yaml file.

```
conda env create -f env/nextflow.yaml
```

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

Edit the `nextflow.config` file to contain the correct output directories and point to the sample sheet in `sample_sheet`, and the references in `fasta_file`, `gtf`, `star_genome_lib`, `cicero_genome_lib`. 

Then  call the wrapper script on the command line, `main_run.sh`. The nextflow config file can be specified by upating the variable `NFX_CONFIG` for the workflow to be run locally, using Singularity/Apptainer on an HPC system. 

To select different executors and container engines, use of the nextflow profiles ('local_apptainer', 'PBS_apptainer', 'local_singularity', 'PBS_singularity'). Update the `NFX_PROFILE` in main_run.sh to select a different executor/container engine. 

```
conda activate nextflow
./main_run.sh
```

# CICERO Fusion Detection

This repo also provides a NF workflow for the [CICERO](https://github.com/stjude/CICERO) fusion detection algorithm. The genome references for CICERO can be found [here](https://www.bcgsc.ca/downloads/genomes/9606/hg19/1000genomes/bwa_ind/genome/GRCh37-lite.fa) and [here](https://doi.org/10.5281/zenodo.3817656).

For [GRCh38 fasta](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)


## Quality Control 

In addition, running [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [multiQC](https://multiqc.info/) on the input fastq files. 

# DSL2 Resources

* [Tutorial](https://github.com/nextflow-io/nfcamp-tutorial)
* [Documentation](https://www.nextflow.io/docs/latest/dsl2.html)
* [Dynamic memory allocation](https://lucacozzuto.medium.com/handling-failing-jobs-with-nextflow-24405b97b679)

# Resources on Automated Builds
* [Fred Hutch SciWiki](https://sciwiki.fredhutch.org/hdc/hdc_building_containers/#step-2-define-your-docker-image)