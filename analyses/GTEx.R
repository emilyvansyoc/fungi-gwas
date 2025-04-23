## updated GTEx with FAVs
# 12/2024
### 1/2/2025 update with GTEx v10 

library(tidyverse)
library(arrow) # to read in v10 parquet files

# get FAVs
ids <- read.table("data/unique_FAVs.txt", header = TRUE)
load("data/Nov2024_sigs_SNPandStructural.RData")
sigs <- newsigs %>% mutate(chr = paste0("chr", X.CHROM))

## ----- convert rsIDs into b38 genomic coordinates -----

## get GTEx lookup table (too big to upload to Github)
# available from GTEx portal 
gtex <- read.table("GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt", sep = "\t", header = TRUE)

# decrease size 
g1 <- gtex %>% filter(chr %in% sigs$chr)
rm(gtex)

# format sigsnps for matching
sigs <- sigs %>% 
  mutate(variant_id_b37 = paste(X.CHROM, POS, REF, ALT, "b37", sep = "_")) %>% 
  mutate(pos.flip = paste(X.CHROM, POS, ALT, REF, "b37", sep = "_"))

# match by SNPID
m1 <- g1 %>% 
  inner_join(sigs %>% dplyr::select(ID, pos.flip, myvariant_id_b37 = variant_id_b37), by = c("rs_id_dbSNP151_GRCh38p7" = "ID"))
# 145 of 148 unique variants 

nom1 <- sigs %>% dplyr::select(ID, variant_id_b37, pos.flip) %>% anti_join(g1, by = c("ID" = "rs_id_dbSNP151_GRCh38p7")) # four; two are structural (not in GTEx)
# these don't exist in GTEx - confirmed with manual lookup in the database (including swapped positions of ref and alt)

## from the matches, verify that the ref/alt alleles are the same and get b38 IDs
m1 <- m1 %>% 
  mutate(matches = case_when(
    variant_id_b37 == myvariant_id_b37 ~ "perfect",
    variant_id_b37 == pos.flip ~ "flipped"#,
    #variant_id_b37 != myvariant_id_b37 & variant_id_b37 != pos.flip ~ "nomatch"
  )) 
unique(m1$matches) # all are either a perfect match or have flipped alleles
length(unique(m1$variant_id)) # 144

# get "flipped" positions to match gtex
m1 <- m1 %>% 
  mutate(myvariant_id_b37_correctorder = if_else(matches == "perfect", myvariant_id_b37, pos.flip))

#### ---- get GTEx data ----

## search for presence of these SNPs in the GTEx data
# downloaded from: https://gtexportal.org/home/downloads/adult-gtex/qtl
### 1/2/2025 update - use recently released v10
### NOT shared on Github since this is publicly available and big

## get list of significant eQTL-variant pairs
outdf1 <- data.frame()
mypath <- "path/GTEx_v10/GTEx_Analysis_v10_eQTL_updated"
myfi <- list.files(mypath, pattern = "*signif*", full.names = TRUE)
for(i in 1:length(myfi)) {
  cat("working on", i, "of", length(myfi), "\n")
  fi <- read_parquet(myfi[i])
  fi$tissue <- str_remove(myfi[i], "path/GTEx_v10/GTEx_Analysis_v10_eQTL_updated/")
  fi$tissue <- str_remove(fi$tissue, ".v10.eQTLs.signif_pairs.parquet")
  outdf1 <- rbind(outdf1, fi)
}
# save 
save(outdf1, file = "path/Nov2024_GTEx_sigSNPs_fromLUtable_v10.RData")


#### ---- get eQTLs for significant snps ----


eqtl <- outdf1 %>% filter(variant_id %in% m1$variant_id) 
# IN V10; 102 (82 unique variants) - 5 unique tissues including liver 
eqtl2 <- outdf1 %>% filter(variant_id %in% m1$pos.flip) # none with the flipped position

# join to our data
eqtl1 <- eqtl %>% full_join(m1, relationship = "many-to-many") 

length(unique(eqtl$variant_id[!is.na(eqtl$gene_id)])) # 82 unique variants

# save
save(eqtl, file = "data/Nov2024_GTEx_sigSNPs_fromLUtable_v10.RData")
save(eqtl1, file = "data/Nov2024_GTEx_allFAVs_v10.RData")

## ---- investigate eQTLs ----

unique(eqtl$variant_id)
unique(eqtl$gene_id) # 5 genes; in v10, now 7 genes
unique(eqtl$tissue) # three tissues; in v10, now 6 tissues

### which taxa are these associated with?
load("data/Nov2024_sigs_SNPandStructural.RData")

# get sigs and fungal taxa
sigIDs <- unique(eqtl1$rs_id_dbSNP151_GRCh38p7[!is.na(eqtl1$gene_id)])
sigs <- newsigs %>% filter(ID %in% sigIDs) %>% dplyr::select(taxa, ID)
length(unique(sigs$taxa)) # 5 taxa; Candida, Sacch family, Aspergillaceae, Kazachstania, and Capnodiales
sigs %>% group_by(taxa) %>% count() # Kaz has 35, Capnodiales has 45

## to keep things straight, add the two fungal taxa with the same SNPID together 
sigs1 <- sigs %>% 
  mutate(taxa = if_else(str_detect(taxa, "Candida"), "Genus_Candida_AND_Family_Saccharomycetales_fam_Incertae_sedis", taxa)) %>% 
  filter(!taxa == "Family_Saccharomycetales_fam_Incertae_sedis")

# join together and investigate
full <- eqtl1 %>% drop_na() %>% full_join(sigs1, by = c("rs_id_dbSNP151_GRCh38p7" = "ID")) %>% distinct()

# add gene name by hand (look up in Ensemble)
full <- full %>% 
  mutate(gene_name = case_when(
    gene_id == "ENSG00000140945.17" ~ "CDH13",
    gene_id == "ENSG00000155511.18" ~ "GRIA1",
    gene_id == "ENSG00000164164.17" ~ "OTUD4",
    gene_id == "ENSG00000166833.23" ~ "NAV2",
    gene_id == "ENSG00000164161.10" ~ "HHIP",
    gene_id == "ENSG00000248890.3" ~ "HHIP-AS1",
    gene_id == "ENSG00000250539.1" ~ "KRT8P33"
  ))

full %>% group_by(taxa, gene_name, tissue) %>% count()
# all 45 of the FAVs in Capnodiales are in the same gene
# Kaz remains the only taxa with eQTLs in multiple genes/tissues (4 different genes)

# save
save(full, file = "data/Nov2024_GTEx_fullFAVResults_v10.RData")
# save for supplementary table
forsupp <- full %>% dplyr::select(gene_id, gene_name, rs_id_dbSNP151_GRCh38p7, tissue, variant_id, ref, alt, tss_distance, af, pval_nominal, slope)
write.table(forsupp, file = "data/results_GTEx_forSupplementary.txt", sep = "\t", row.names = FALSE)

# how many variants are associated with more than one gene?
multis <- full %>% group_by(variant_id) %>% count() %>% filter(n>1) # 10 are associated with all 3 genes; all in same genomic region
full %>% filter(variant_id %in% multis$variant_id) %>% group_by(taxa, gene_name) %>% count() # all are Kaz, HHIP, HHIP-AS1, and OTUD4

## ---- which eQTLs are also overlapping the gene of interest? ----

# get SNPNexus results
genes <- read.table("data/SNPNexus_results/near_gens_allSNPNov2024.txt", sep = "\t", header = TRUE)

# filter and join
overlap <- full %>% inner_join(genes, by = c("rs_id_dbSNP151_GRCh38p7" = "Variation.ID", "gene_name" = "Overlapped.Gene")) %>% 
  filter(!gene_name == "None")
# which are these?
length(unique(overlap$variant_id))
overlap %>% group_by(taxa, gene_name) %>% count()
