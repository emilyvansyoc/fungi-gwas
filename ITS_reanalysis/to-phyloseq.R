## get OTU table and sintax information and output phyloseq object
# EVS 12/2024

suppressMessages(library(phyloseq))
suppressMessages(library(microViz))
suppressMessages(library(vegan))
suppressMessages(library(EcolUtils))
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

## set arguments 
args <- commandArgs(TRUE)
otutab <- args[[1]] # OTU file from VSEARCH output
taxtab <- args[[2]] # SINTAX file from VSEARCH output
refseqs <- args[[3]] # OTU sequence in fasta file
outpath <- args[[4]] # pathway to save output files
acc <- args[[5]] # accession (for naming)

### ---- format for phyloseq ----

# OTU table
otab <- read.table(otutab, sep = "\t", header = TRUE, comment.char = "", na.strings = "") %>% column_to_rownames(var = "X.OTU.ID")

# SINTAX (wonky)
sin <- read.table(taxtab, sep = "\t", header = FALSE, comment.char = "", na.strings = "") %>%  dplyr::select(V1, V4) %>% separate_wider_delim(cols = V4, names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), delim = ",", too_few = "align_start") %>% mutate(OTU_ID = str_extract(V1, "OTU_(\\d){1,10}")) %>% mutate(across(!OTU_ID, ~str_remove(.x, "[:alpha:]:"))) %>% column_to_rownames(var = "OTU_ID") %>% dplyr::select(-V1)

# reference sequences
seqs <- readDNAStringSet(filepath = refseqs)
names(seqs) <- str_extract(names(seqs), "OTU_(\\d){1,1000}")

# make phyloseq
ps <- phyloseq(otu_table(otab, taxa_are_rows = TRUE),
               tax_table(sin %>% as.matrix()))
ps@refseq <- seqs

## save raw phyloseq object
save(ps, file = paste0(outpath, "/raw_phyloseq_", acc, ".RData"))

## ---- rarefy and filter ----

# remove non-Fungi and NA phyla
ps1 <- ps %>% tax_select("Fungi", "Kingdom")
ps2 <- ps1 %>% tax_fix() %>% tax_select("Fungi Kingdom", "Phylum", deselect = TRUE)

## rarefy
otab <- ps2 %>% otu_get() %>% as.data.frame()

# rarefy at 1500
otab1 <- otab[rowSums(otab) > 1500, ]
rtab <- suppressMessages(EcolUtils::rrarefy.perm(otab1, sample = 1500, n = 1000, round.out = TRUE))

# re-phyloseq
rarephy <- phyloseq(
  otu_table(rtab, taxa_are_rows = FALSE),
  tax_table(ps2)
)

# remove empty taxa
rarephy <- rarephy %>% tax_filter(min_prevalence = 1, min_total_abundance = 1)

# add reference sequences
rarephy@refseq <- ps2@refseq

## ---- save ----
psf <- rarephy
save(psf, file = paste0(outpath, "/rarefied_phyloseq_", acc, ".RData"))
