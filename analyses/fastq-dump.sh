#!/bin/bash

# run locally; get reads from SRA

# change output variables in terminal, not in script

### activate micromamba
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate bioinfo
set -uex

### ---- set variables ----

# directory to run ID's, include name of runid file
RUNIDS=$1

# directory to output files
OUTDIR=$2

# directory for quality analysis
QUALDIR=$3


## ---- get reads from SRA with fastq-dump ---

# make output directory if it doesn't exist
mkdir -p $OUTDIR

# fastq-dump
cat $RUNIDS | head -50 | parallel /Users/epb5360/bin/sratoolkit.3.0.0-mac64/bin/fastq-dump {} --outdir $OUTDIR

echo "read dump complete"

## ---- run fastqc and multiqc for quality ----

# make output directory if it doesn't exist
mkdir -p $QUALDIR

# run fastqc on forward reads 
cat $RUNIDS | head -50 | parallel fastqc $OUTDIR/{}_1.fastq -o $QUALDIR -q 

# run fastqc on reverse reads 
#cat $RUNIDS | head -50 | parallel fastqc $OUTDIR/{}_2.fastq -o $QUALDIR

# run multiqc on all reads
multiqc $QUALDIR* -o $QUALDIR -n rawreads.multiqc

echo "quality check complete"

