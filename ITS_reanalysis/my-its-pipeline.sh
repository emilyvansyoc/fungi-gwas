#### master pipeline for pulling ITS data from SRA and re-analyzing 
# EVS 12/2024

source ~/.bashrc
eval "$(micromamba shell hook --shell bash)"
micromamba activate itsenv
set -uex

## ---- CHANGE THESE VARIABLES ----

# accession number
ACC=$1

# input directory
DIR=/path/to/my/working/directory
mkdir -p $DIR/$ACC

# set prefix (European is ERR, American is SRR)
PREF=$2

# location to SINTAX UNITE database
DB=/path/to/my/database

## ---- get runinfo ----

# get runinfo
esearch -db sra -query $ACC | efetch -format runinfo > $DIR/runinfo_$ACC.csv

# get runids
cat $DIR/runinfo_$ACC.csv | cut -f 1 -d ',' | grep $PREF > $DIR/runids_$ACC.txt

# print message
echo "\n Number of samples detected: "
cat $DIR/runids_$ACC.txt | wc -l

## ---- fastq-dump ----

# make directory
RAW=$DIR/$ACC/rawreads
mkdir -p $RAW

# run fastq-dump on all samples
cat $DIR/runids_$ACC.txt | parallel fastq-dump --split-files --outdir $RAW {}

# print stats for the reads
seqkit stat $RAW/*.fastq -T -j 10 -o $DIR/$ACC/stats_raw_$ACC.txt

## ---- run primer removal through dada2 ----

# requires two external arguments: path to raw data, and primer file
# if primer file is not available, this will quit without an error
Rscript dada2-cutadapt.R $RAW $DIR/primers_$ACC.fasta

# now, primer removed files are located in a subdirectory
CUTAD=$RAW/cutadapt

## ---- run OTU generation through VSEARCH ----

### this creates several subdirectories

# location for quality filtered reads
QUAL=$DIR/$ACC/qualfilt/
mkdir -p $QUAL

# merged reads
MERGE=$DIR/$ACC/merge/
mkdir -p $MERGE

# dereplicated reads
DEREP=$DIR/$ACC/derep/
mkdir -p $DEREP

## run the code VSEARCH code ######

## error catch: if primers were not removed with cutadapt, use filtN files instead
testfile=($(ls $RAW | head -1 | cut -f 1 -d _ ))
if [ -e "$CUTAD/${testfile}_1.fastq" ]; then
    echo "VSEARCH MERGING ON CUTADAPTED READS"
    # merge paired-end reads
    cat $DIR/runids_$ACC.txt | parallel vsearch --fastq_mergepairs $CUTAD/{}_1.fastq --threads 10 --reverse $CUTAD/{}_2.fastq --fastq_minovlen 30 --fastq_allowmergestagger --fastqout $MERGE/{}.fastq

    # print stats for the reads
    seqkit stat $CUTAD/*.fastq -T -j 10 -o $DIR/$ACC/stats_cutadapt_$ACC.txt

    else
        echo "NO CUTADAPTED FILES FOUND: VSEARCH MERGING ON FILTN READS"
        # merge paired-end reads
        cat $DIR/runids_$ACC.txt | parallel vsearch --fastq_mergepairs $RAW/filtN/{}_1.fastq --threads 10 --reverse $RAW/filtN/{}_2.fastq --fastq_minovlen 30 --fastq_allowmergestagger --fastqout $MERGE/{}.fastq
fi

# run stats
seqkit stat $MERGE/*.fastq -T -j 10 -o $DIR/$ACC/stats_merged_$ACC.txt

# quality filtered merged reads
cat $DIR/runids_$ACC.txt | parallel vsearch --fastx_filter $MERGE/{}.fastq --fastq_minlen 20 --fastq_truncqual 10 --fastaout $QUAL/{}.fasta

# sample-wise dereplication
cat $DIR/runids_$ACC.txt | parallel vsearch --derep_fulllength $QUAL/{}.fasta --output $DEREP/{}_derep.fasta --sizeout --relabel {}.

# merge all samples
cat $DEREP/*_derep.fasta > $DIR/$ACC/all.fasta 

# dereplicate whole datasets
vsearch --derep_fulllength $DIR/$ACC/all.fasta --minuniquesize 2 --sizein --sizeout --output $DIR/$ACC/all.derep.fasta 

### NOTE: there is NO masking done in clustering or de novo steps; this can be changed

# cluster
vsearch --cluster_size $DIR/$ACC/all.derep.fasta --threads 10 --id 0.97  --qmask none --sizein --sizeout --centroids $DIR/$ACC/centroids.fasta

# sort and remove singletons
vsearch --sortbysize $DIR/$ACC/centroids.fasta --sizein --sizeout --minsize 2 --output $DIR/$ACC/sorted.fasta

# de novo chimeras
vsearch --uchime_denovo $DIR/$ACC/sorted.fasta --sizein --sizeout --qmask none --nonchimeras $DIR/$ACC/denovo.nonchimeras.fasta

# relabel OTUs
vsearch --fastx_filter $DIR/$ACC/denovo.nonchimeras.fasta --sizein --sizeout --relabel OTU_ --fastaout $DIR/$ACC/otus.fasta

# assign OTUs to sequences
vsearch --usearch_global $DIR/$ACC/all.fasta --threads 10 --db $DIR/$ACC/otus.fasta --id 0.97 --strand plus --sizein --sizeout --qmask none --dbmask none --otutabout $DIR/$ACC/otutab.txt

# run SINTAX
vsearch --sintax $DIR/$ACC/otus.fasta --db $DB --tabbedout $DIR/$ACC/sintax50.txt --sintax_cutoff .50

### ---- make into phyloseq object, rarefy and filter ----

# make output path
OUT=$DIR/$ACC/ps-out/
mkdir -p $OUT

# run RScript
Rscript to-phyloseq.R $DIR/$ACC/otutab.txt $DIR/$ACC/sintax50.txt $DIR/$ACC/otus.fasta $OUT $ACC