### separate 'all variants' into SNPs and structural
# EVS 12/2024

library(tidyverse)
library(gridExtra)
library(ggpubr)
#library(qqman)
#library(phyloseq)
#library(microViz)
#library(cowplot)

# set pathway to output files (restricted data)
path <- "restricted/updateITS_results_allsnp_allvariant/"
### all files are read in, concatenated, and saved as the below RData object

### ---- get effect sizes ----

# get effect sizes (restricted data, not shared) - also huge
load("restricted/effectsizes_allSNP_allvariants_alltaxa.RData")

### 11/2024; get FAVs at each significance level (exploratory, genome-wide, study-wide) and combine into one dataset
#### remove Saccharomycetes Class from these tallies because it is identical to the Saccharomycetales Order

# 1. study-wide significance
sigs <- alleffects %>% 
  # remove Saccharomycetes Class
  filter(!str_detect(taxa, "Class_Saccharomycetes")) %>% 
  mutate(studysig = p.adjust(P, method = "fdr"))
# get study-wide significance at Q < 0.2
range(sigs$studysig)
study.sigs <- sigs %>% filter(studysig < 0.2) # 10
unique(study.sigs$X.CHROM) # 2 chromosomes
unique(study.sigs$taxa) # 2 taxa
study.sigs %>% filter(nchar(REF) == 1 & nchar(ALT) == 1) # 9 of 10 are SNPs

# 2. genome-wide significance
gen.sigs <- sigs %>%
  filter(P < 5e-8)
unique(gen.sigs$ID) # 24 unique IDs - one ID is detected for two fungi (rs6879769)
unique(gen.sigs$X.CHROM) # 4 chromosomes
unique(gen.sigs$taxa) # 7 taxa
gen.sigs %>% filter(nchar(REF) == 1 & nchar(ALT) == 1) %>% nrow() # 23 of 24 are SNPs

# 3. exploratory significance (FDR correction for variants only)
sigs <- sigs %>% 
  group_by(taxa) %>% 
  mutate(expsig = p.adjust(P, method = "fdr")) %>% 
  ungroup() 
ex.sigs <- sigs %>% 
  filter(expsig < 0.05) # 143 unique IDs
# how many of these are "only" exploratory (i.e. don't meet genwide or study wide sig?)
ex.sigs %>% 
  mutate(genwide = if_else(P < 5e-8, TRUE, FALSE),
         issnp = if_else(nchar(REF) == 1 & nchar(ALT) == 1, TRUE, FALSE)) %>% 
  group_by(genwide, issnp) %>% count() # 124 "additional" sigs 
unique(ex.sigs$X.CHROM)
unique(ex.sigs$taxa)


## there are 6 SNPs that DO meet genome-wide sig that do NOT meet the exploratory FDR; this makes the math wonky
#### 4. get ALL sigs 
allsigs <- sigs %>% 
  filter(P < 5e-8)  %>% 
  rbind(ex.sigs) %>% distinct()
length(unique(allsigs$ID)) # 148 unique IDs
unique(allsigs$X.CHROM) # 7 chromosomes
unique(allsigs$taxa) # 9 taxa 
allsigs %>% filter(nchar(REF) == 1 & nchar(ALT) == 1) %>% nrow() # 140 are SNPs

## format
newsigs <- allsigs %>% 
  dplyr::select(taxa, ID, X.CHROM, POS, BETA, A1_FREQ, REF, ALT, P, expsig, studysig) %>% 
  mutate(taxa = str_remove(taxa, "allsnp_allvariants\\."),
         taxa = str_remove(taxa, "\\.glm\\.linear")) %>% 
  mutate(posID = paste(X.CHROM, POS, sep = ":")) %>% 
  #dplyr::select(-c(X.CHROM, POS)) %>% 
  relocate(taxa, ID, posID) %>% 
  arrange(P) %>% ungroup()

# save
save(newsigs, file = "data/Nov2024_sigs_SNPandStructural.RData")

# get just the FAVs
fav <- newsigs %>% ungroup() %>% dplyr::select(ID) %>% distinct()
write.table(fav, file = "data/unique_FAVs.txt", sep = "\t", row.names = FALSE)

# format for SNPNexus
nex <- fav %>% mutate(Type = "dbsnp") %>% rename(Name = ID) %>% relocate(Type)
write.table(nex, file = "data/unique_FAVs_forSNPNexus.txt", sep = "\t", row.names = FALSE, quote = FALSE)
