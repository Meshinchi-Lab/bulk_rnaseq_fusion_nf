# STAR-fusion NF

 STAR-fusion run on AWS Batch using the Nextflow workflow manager.

 This workflow requires:
  1. Sample sheet, tab delimited
  2. The location of the CTAT resource library

The output of STAR aligner, STAR-Fusion, and Fusion-Inspector will be uploaded to an S3 bucket. This includes the most relevant output files, such as SJ.out.tab, aligned.bam, and chimeric.junctions.tab, and the fusion inspector HTML report.

The output files will be put into a directory that is named after the sample ID provided in the sample sheet file.  

# To Run

## First Step:

First, create a sample manifest for the fastq files that are hosted in an S3 bucket. The manifest file is a simple 3 column, tab-delimited text file with the column names "Sample", "R1", and "R2", where R1 and R2 are paired end-fastqs.

I have created a sample manifest file that can be used to select the appropriate files, with a demonstration in `Nextflow_AWS_Sample_Sheets_from_Manifest.ipynb` that uses the associated `create_sample_sheet.py`. Below is a very simple example directly on the command line.

```
python3 create_sample_sheet.py "fh-pi-meshinchi-s" "fastq" --prefix "SR/picard_fq2/"
```

## Second Step

Edit the `STAR_Fusion_run.sh` file to contain the correct output directory in `--output_folder` and point to the correct sample sheet in `--sample_sheet`. You can also update the filename for the output html report in `-with-report`. Then just simple call the script on the command line.

note: I have access to an HPC with specific software modules that can be loaded. If you have a custom installation of Nextflow, simply make sure the nextflow executable is in your PATH.

```
./STAR_Fusion_run.sh
```

# CICERO Fusion NF

This repo also provides NF workflow for the [CICERO](https://github.com/stjude/CICERO) fusion detection algorithm. The image was build from the Dockerfile provided on the CICERO github without modification for v.0.3.0 and v.0.2.0.  

The next steps are to create a [DSL2 Nextflow](https://www.nextflow.io/docs/latest/dsl2.html) workflow in order to make this a multipart workflow with options for runnng STAR or CICERO or both.
