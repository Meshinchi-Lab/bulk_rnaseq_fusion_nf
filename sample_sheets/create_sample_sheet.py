#!/usr/bin/env python3

#Jenny Smith
#June 24, 2020 
#Purpose: create a sample sheet of files hosted on AWS S3. 

#import modules 
import argparse
import os
import re
import itertools
import boto3
import logging
from botocore.exceptions import ClientError
import numpy as np
import pandas as pd


#set-up command line arguments
parser = argparse.ArgumentParser()

parser.add_argument('bucket',
                    help='''AWS S3 bucket name. 
                    No slashes or S3:// included.
                    Example fh-pi-my-bucket''')

parser.add_argument('filetype',default="fastq",
                    help='''The file type as a string.
                    Accepted inputs are either fastq or bam. 
                    The default is fastq.''')

parser.add_argument('--filename', default="sample_sheet.txt",
                    help='''A string of sample IDs that
                    are part of the filenames being processed.
                    Can be comma separated or space separated.''')

parser.add_argument('--prefix', default="",
                    help='''The AWS S3 prefix to search for the bam or
                    fastq files. The delimiter is assumed to be 
                    '/' forward slash. Must have a trailing slash.''')

parser.add_argument('--samples', default="",
                    help='''A string of sample IDs that
                    are part of the filenames being processed.
                    Can be comma separated or space separated.''')

args = parser.parse_args()



#Connection to S3 Bucket 
s3 = boto3.resource('s3')
bucket = s3.Bucket(args.bucket)
prefix = args.prefix
delim = "/"



#Processes for paired end fastqs 
if args.filetype == "fastq":
   
    #function to parse the object summary from Boto3 for fastqs
    def sample_name(s3_object_summary):
        sample = s3_object_summary.key.split("/")[2]
        sample = re.sub(r"_[Rr][12].+$","", sample)
        return(sample)

    #query S3 bucket to list the fastq files 
    fqs = bucket.objects.filter(Delimiter=delim,
                                Prefix=prefix)
    


    #iterate over the fastqs 
    fqs_PE = dict()
    for obj in fqs:
        samp = sample_name(obj)
        
        #Some of the listed files are just the prefix, 
        #and some are stderr files for processing. 
        #So skip these types of files 
        if re.search(r".stderr$",samp) or len(samp)==0:
            continue

        if re.search(r"[_-][Rr]1", obj.key):
            r1 = '{0}//{1}/{2}'.format("s3:", obj.bucket_name, obj.key)

            if samp not in fqs_PE.keys():
                fqs_PE[samp] = [r1]
            else:
                fqs_PE[samp] = fqs_PE[samp] + [r1]

        elif re.search(r"[_-][Rr]2", obj.key):
            r2 = '{0}//{1}/{2}'.format("s3:", obj.bucket_name, obj.key)

            if samp not in fqs_PE.keys():
                fqs_PE[samp] = [r2]
            else:
                fqs_PE[samp] = fqs_PE[samp] + [r2]
        
        #finally sort just to ensure correct order
        fqs_PE[samp].sort()

    #Create a pandas dataframe and save it 
    #split the samples for filtering by space or comma
    string = args.samples
    fname = args.filename
    
    regex = re.compile('|'.join(re.split(',| ',string)))
    filtered = [{"Sample":sample,"R1":fastqs[0],"R2":fastqs[1]} for sample, fastqs in fqs_PE.items() if re.search(regex, sample)]
    
    sample_sheet = pd.DataFrame(filtered) 
    sample_sheet.to_csv(path_or_buf=fname,
                    sep="\t", header=True,
                    index=False,quoting=None)