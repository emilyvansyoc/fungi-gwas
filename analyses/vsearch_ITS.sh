#!/bin/bash

### run VSEARCH locally under a micromamba install
# EVS 10/2023, updated 1/2024 remove quality filter (do in R)

# activate environment
source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate vsearch
set -uex

# location to raw reads
RAW=my-path/filtered

# working directory
DIR=my-path/

# set location to file IDs
IDS=$DIR/runids_HMP.txt

# set location to put merged reads  
MERGE=$DIR/merged/
mkdir -p $MERGE

# set location to put merged-to-fasta reads
TOFA=$DIR/merged-to-fasta
mkdir -p $TOFA

# set location to put dereplicated reads
DEREP=$DIR/derep/
mkdir -p $DEREP


### ----- run code -----

# merge paired end reads
cat $IDS | parallel vsearch --fastq_mergepairs $RAW/{}_R1_001-ITS.fastq.gz --threads 2 --reverse $RAW/{}_R2_001-ITS.fastq.gz --fastq_minovlen 10 --fastq_maxdiffs 20 --fastq_allowmergestagger --fastqout $MERGE/{}_merge.fastq --fastq_eeout --eetabbedout $DIR/stats_vmerge.txt

# convert to fasta
cat $IDS | parallel seqkit fq2fa $MERGE/{}_merge.fastq -o $TOFA/{}.fasta

# do sample-wise dereplication
cat $IDS | parallel vsearch --derep_fulllength $TOFA/{}.fasta --output $DEREP/{}_derep.fasta --sizeout --relabel {}.

# merge all samples
cat $DEREP/*_derep.fasta > $DIR/all.fasta

# dereplicate whole dataset
vsearch --derep_fulllength all.fasta --minuniquesize 2 --sizein --sizeout --output all.derep.fasta --uc all.derep.uc

## ----- cluster at 97% and then remove chimeras & singletons ----

### NOTE: there is NO masking done in clustering or de novo steps; this can be changed

# cluster
vsearch --cluster_size all.derep.fasta --threads 10 --id 0.97 --strand plus --qmask none --sizein --sizeout --centroids centroids.fasta

# sort and remove singletons
vsearch --sortbysize centroids.fasta --sizein --sizeout --minsize 2 --output sorted.fasta

# de novo chimeras
vsearch --uchime_denovo sorted.fasta --sizein --sizeout --qmask none --nonchimeras denovo.nonchimeras.fasta

# relabel OTUs
vsearch --fastx_filter denovo.nonchimeras.fasta --sizein --sizeout --relabel OTU_ --fastaout otus.fasta

# assign OTUs to sequences
vsearch --usearch_global all.fasta --threads 10 --db otus.fasta --id 0.97 --strand plus --sizein --sizeout --qmask none --dbmask none --otutabout otutab.txt