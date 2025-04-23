## Pheweb results
# EVS 3/2024

library(tidyverse)
library(ggpubr)
library(ieugwasr)
library("TwoSampleMR")

## ---- parse pheWEB results ----

### each FAV was searched by hand in the PheWeb API 
# get results from text file
dat <- read.table("data/pheweb_results.txt", sep = "\t", header = TRUE) %>% 
  dplyr::select(SNPID, Category, Phenotype, P_value, effect_size, SE, Number_of_samples) %>% 
  drop_na()
## filter to FAVs
myids <- read.table("data/unique_FAVs.txt", sep = "\t", header = TRUE) 
dat <- dat %>% 
  filter(SNPID %in% myids$SNPID) # 24 of the original 25 


# get the phenotypes that are above the 7.05e-5 cutoff
### Bonferonni correction for 1419 total codes in the Phewas database
bonf05 <- 0.05/1419
bonf1 <- 0.1/1419
# get strictest sigs
sigs <- dat %>% filter(P_value < bonf05) # 14

# get taxa and match to SNPs
# get list of taxa
load("data/Nov2024_sigs_SNPandStructural.RData")
tax <- allsig %>% filter(!Taxa %in% "Class_Saccharomycetes") %>% dplyr::select(Taxa, SNPID, `Variant Type`) %>% drop_na()

## join
sigs <- sigs %>% left_join(tax)

# summary stats
summary(sigs$effect_size)
summary(sigs$SE)

### ---- validate with ieugwasr package ----

### website to get authentication token (login through Github):
# https://api.opengwas.io/profile/

## set a character vector to the authentication token
mytoken <- "mytoken"


# get SNP list
myids <- read.table("data/unique_FAVs.txt", sep = "\t", header = TRUE) 
snplist <- myids$SNPID

# perform phewas
# get associations from a particular study (fast)
associations(variants = c("rs12149890", "rs12929586"),
             id = c("bbj-a-159", "ebi-a-GCST003116"))

# what does 1000 Genomes say about these variants?
ieugwasr::afl2_rsid(rsid = c("rs12149890", "rs12929586"), opengwas_jwt = mytoken)


