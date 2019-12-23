# STAR-fusion NF

 STAR-fusion run on AWS Batch using the Nextflow workflow manager.

 This workflow requires:
  1. Sample sheet, tab delimited
  2. The location of the CTAT resource library

The output of STAR aligner, STAR-Fusion, and Fusion-Inspector will be uploaded to an S3 bucket. This includes the most relevant output files, such as SJ.out.tab, aligned.bam, and chimeric.junctions.tab, and the fusion inspector HTML report. 

The output files will be put into a directory that is named after the sample ID provided in the sample sheet file.  
