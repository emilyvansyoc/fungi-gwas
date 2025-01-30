#!/bin/bash

### adapted from the Biostar handbook
## run locally

## script to get runinfo and test accession number

# activate conda
conda activate bioinfo

# catch errors
set -uex

## ---- CHANGE THESE VARIABLES ----

# accession number (starts with PRJNA)
ACC=$1

# input directory
DIR=$2

# general output name
OUTNAME=$3

# set prefix (European is ERR, American is SRR)
PREF=SRR

## ---- workflow: get runinfo and test accession ----

# get runinfo
esearch -db sra -query $ACC | efetch -format runinfo > $DIR/runinfo_$OUTNAME.csv

# get runids
cat $DIR/runinfo_$OUTNAME.csv | cut -f 1 -d ',' | grep $PREF > $DIR/runids_$OUTNAME.txt

# get one sample with 10,000 reads to test connection
cat $DIR/runids_$OUTNAME.txt | head -1 | parallel /storage/work/epb5360/miniconda3/envs/bioinfo/bin/fastq-dump -X 10000 --split-files --outdir $DIR {}

# print stats for the reads
seqkit stat $DIR/*.fastq

# print finished message
echo "done!"

