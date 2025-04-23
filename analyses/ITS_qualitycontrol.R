# ITS preliminary quality control
# EVS 8/2023, updated 1/2024 to match hominid methods

library(dada2)
library(tidyverse)
library(ShortRead)
library(Biostrings)
library(phyloseq)

# set pathway - this is a dummy pathway (raw files not provided in Github)
path <- "path-to-raw-files/"
list.files(path)
# how many total?
length(list.files(path, pattern = "_1.fastq")) # 390 total

## ---- remove 18S samples ----

# get list of 18S samples to remove 
# this is the default SRA metadata file
its <- read.table("mypath/SRAMetadata_Nash2017_PRJNA356769.txt", sep = "\t", header = TRUE) %>% 
  filter(str_detect(`Library Name`, "ITS")) %>% 
  # jankify
  mutate(full = paste0("/path//", Run, "_1.fastq")) %>% 
  mutate(fullrev = paste0("path//", Run, "_2.fastq"))

# stupid dada2 part
fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnFs <- fnFs[fnFs %in% its$full]
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))
fnRs <- fnRs[fnRs %in% its$fullrev]

### ---- remove primers ----- 

### get primer sequences from the paper
# ITS3
FWD = "GCATCGATGAAGAACGCAGC"
# ITS4
REV = "TCCTCCGCTTATTGATATGC"

# make sure we get the correct orientation of the primers

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  #require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, 
               Complement = Biostrings::complement(dna), 
               Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

## remove N's from the samples to make primer searching easier
#fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
#fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
#filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, verbose = TRUE)

# once that has been run once, set path file for N-removed files
npath <- "path/filtN/"
fnFs.filtN <- sort(list.files(npath, pattern = "_1.fastq", full.names = TRUE))
fnRs.filtN <- sort(list.files(npath, pattern = "_2.fastq", full.names = TRUE))


# get number of primer hits
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# get primers for all reads
for(i in 1:length(fnFs.filtN)) {
  
  # print sample name
  cat(paste0("\n sample ", fnFs.filtN[i], "\n"))
  
  print(
    # print primers
    rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[i]]), 
          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[i]]),
          REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[i]]), 
          REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[i]]))
  )
  
}

# from this;  have primers in the revcomp orientation in FWD.REverseReads and REV.ForwardReads
# anywhere from a few to several thousand 
# one sample has Forward with REV.REverseReads (two samples)
# one has Forward in REV.ForwardReads and REV.ReverseReads

# so, to remove these;
# will search for the forward primer in the forward reads in normal orientation
# search for reverse primer in reverse complement in the forward and reverse reads


## ---- remove primers with cutadapt in R not shell (following tutorial) ----
cutadapt <- "path-to/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs.filtN))
fnRs.cut <- file.path(path.cut, basename(fnRs.filtN))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs.filtN)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

## sanity check; ensure primers have been removed
for(i in 1:length(fnFs.cut)) {
  
  # print sample name
  cat(paste0("\n sample ", fnFs.cut[i], "\n"))
  
  print(
    # print primers
    rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[i]]), 
          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[i]]),
          REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[i]]), 
          REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[i]]))
  )
  
}

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)

# plot quality
plotQualityProfile(fnRs.cut[1])



# ---- filter and trim ----

filtFs <- file.path(path, "filtered_Phred30", basename(fnFs.cut))
filtRs <- file.path(path, "filtered_Phred30", basename(fnRs.cut))

# filter and trim
# max N=0
#truncQ=2
#rm.phix=T
# maxEE=2
trimout <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs,
                         truncQ = 30,
                         multithread = TRUE, verbose = TRUE)

## remove empty files
empties <- rownames(trimout)[trimout[,2] == 0]
empties <- paste0("path/filtered_Phred30/", empties)

noem <- filtFs[!filtFs %in% empties] # 3 total
emptyR <- gsub("_1", "_2", empties)
noemR <- filtRs[!filtRs %in% emptyR]


