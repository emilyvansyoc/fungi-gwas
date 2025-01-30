## vsearch-to-phyloseq and prep for GWAS
# EVS 2/2024

library(tidyverse)
library(phyloseq)
library(ggpubr)
library(vegan)
library(EcolUtils)
library(microViz)
library(Biostrings)
library(car)

## ---- build phyloseq object ----

# OTU table
otab <- read.table("data/hmpITS_otutab.txt", sep = "\t", header = TRUE, comment.char = "", na.strings = "") %>% column_to_rownames(var = "X.OTU.ID")

# SINTAX taxonomy
sin <- read.table("data/hmpITS_sintax50.txt", sep = "\t", header = FALSE, comment.char = "", na.strings = "") %>%  dplyr::select(V1, V4) %>% separate_wider_delim(cols = V4, names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), delim = ",", too_few = "align_start") %>% mutate(OTU_ID = str_extract(V1, "OTU_(\\d){1,4}")) %>% mutate(across(!OTU_ID, ~str_remove(.x, "[:alpha:]:"))) %>% column_to_rownames(var = "OTU_ID") %>% dplyr::select(-V1)

# metadata
ids <- readxl::read_xlsx("data/SRAMetadata_Nash2017_PRJNA356769.xlsx") %>% 
  filter(!str_detect(`Library Name`, "18S")) %>% 
  select(Run, `Sample Name`) %>% 
  mutate(id = sapply(str_split(`Sample Name`, "_"), `[`, 1)) %>% 
  mutate(visit = sapply(str_split(`Sample Name`, "_"), `[`, 2)) %>% 
  column_to_rownames(var = "Run") %>% 
  mutate(idvisit = paste(id, visit, sep = "_"))

# reference sequences
seqs <- readDNAStringSet(filepath = "data/otus.fasta")
names(seqs) <- str_extract(names(seqs), "OTU_(\\d){1,1000}")

## make phyloseq
ps <- phyloseq(otu_table(otab, taxa_are_rows = TRUE),
               tax_table(sin %>% as.matrix()),
               sample_data(ids))
ps@refseq <- seqs

# save
save(ps, file = "data/raw-phylo-hmpITS.RData")

## ---- filter fungal kingdom ----

# get only Fungal kingdom
ps1 <- ps %>% tax_select("Fungi", "Kingdom") # 622

# get NA's
nas <- ps1 %>% tax_fix() %>% tax_select("Fungi Kingdom", "Phylum") # 42
taxa_names(nas)
# resolve these after filtering

## ---- filter GWAS ----

##### identical to older methods ###

# get IDs of fungal data
fun_ids <- unique(ps1@sam_data$id)

# get IDs of matched WGS data (has already been matched to fungal IDs)
fam <- read.table("private/final_indepSNP_remhet_sub.fam", stringsAsFactors = FALSE) %>% 
  dplyr::select(V1, V2)
names(fam) <- c("FID", "IID")

# get the IDs of genome data
ids <- read.table("private/hmpIDs.txt", sep = "\t", header = FALSE)
names(ids) <- "GID"

# these are SAMPLE NAMES; get key to match them to RANSID (what the fungi "id" column is)
key <- read.table("private/HMP_samplekey.txt", sep = "\t", header = TRUE)
length(which(key$SN %in% ids$GID))
key <- key %>% filter(SAMPLE_USE == "Seq_DNA_WholeGenome; Seq_DNA_SNP")
# get just the columns we need and samples that match with ITS
keysub <- key %>% 
  dplyr::select(SN, RANDSID) %>% 
  dplyr::rename(IID = SN) %>% 
  inner_join(fam) %>% 
  mutate(id = as.character(RANDSID))

## subset ITS data
pssub <- ps1 %>% 
  ps_filter(id %in% keysub$RANDSID) %>% 
  # for most individuals; this is visit 1
  ps_filter(visit == 1)
# add the two individuals who don't have visit 1 but who have visits 2 or 3
#sub1 <- ps1 %>% 
#ps_filter(id %in% c("788731474", "246515023"))
#pssub <- merge_phyloseq(pssub, sub1)


# rename samples
pssub <- pssub %>% ps_join(keysub)
sample_names(pssub) <- pssub@sam_data$IID


### ---- RAREFY -----

# get data
otab <- pssub %>% otu_get() %>% as.data.frame()
summary(rowSums(otab))

# plot
rarecurve(otab, step = 1000)

# rarefy at different levels
otab1 <- otab[rowSums(otab) > 1500, ]
rtab <- EcolUtils::rrarefy.perm(otab1, sample = 1500, n = 1000, round.out = TRUE)

# plot
rarecurve(rtab, step = 100)

# re-phyloseq
rarephy <- phyloseq(
  otu_table(rtab, taxa_are_rows = FALSE),
  tax_table(pssub),
  sample_data(pssub)
)

# remove empty taxa
rarephy <- rarephy %>% tax_filter(min_prevalence = 1, min_total_abundance = 1)
# remove empty samples
rarephy <- rarephy %>% ps_filter(sample_sums(rarephy) > 5)

## ---- resolve NA taxa before collapsing ----


# get NAs
nas <- rarephy %>% tax_fix() %>% tax_select("Fungi Kingdom", "Phylum") # 30
taxa_names(nas)
# 1206: conflict past order - k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales
# 183: conflict past class - k__Fungi;p__Ascomycota;c__Leotiomycetes
# 206: conflict past class - k__Fungi;p__Ascomycota;c__Sordariomycetes
# 225: k__Fungi;p__Basidiomycota;c__Malasseziomycetes;o__Malasseziales;f__Malasseziaceae;g__Malassezia
# 239: conflict past fam -k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Cordycipitaceae
# 261: UNITE only - k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales
# 281: k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Phaffomycetaceae;g__Starmera
# 3: k__Fungi;p__Basidiomycota;c__Malasseziomycetes;o__Malasseziales;f__Malasseziaceae;g__Malassezia
# 300: k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Auriculariales;f__Exidiaceae
# 350: k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus
# 373: k__Fungi;p__Ascomycota
# 375: k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Trichosphaeriales;f__Trichosphaeriaceae;g__Nigrospora
# 389: k__Fungi;p__Ascomycota
# 395: k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Chaetothyriales;f__Cyphellophoraceae;g__Cyphellophora
# 41: BLAST only - k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__Didymellaceae;g__Epicoccum
# 463: too short and none match 100% -> remove
# 479: k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycodaceae
# 486: too short and none match 100% -> remove
# 519: conflicts - k__Fungi;p__Ascomycota
# 553: k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Hymenochaetales;f__Tubulicrinaceae;g__Tubulicrinis
# 655: no good matches; remove
# 696: k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus
# 741: could be Pleosporales Order but conflict with low scores - remove
# 75: lots of conflict past order - k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales
# 758: k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Pichiaceae;g__Pichia
# 782: k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Sordariales;f__Chaetomiaceae;g__Mycothermus
# 8: all BLAST matches to Candida sake but Debaryomycetaceae has higher scores in UNITE - k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales
# 818: k__Fungi;p__Ascomycota
# 940: no good matches and conflict at phylum level; remove
# 943: k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Trichomonascaceae

# get list to remove
torem <- c("OTU_463", "OTU_486", "OTU_655", "OTU_741", "OTU_940")
# get list to keep
torep <- taxa_names(nas)[!taxa_names(nas) %in% torem]

## make dataframe of replacements
newnames <- data.frame(OTU = torep,
                       Tax = c("k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales",
                               "k__Fungi;p__Ascomycota;c__Leotiomycetes",
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes",
                               "k__Fungi;p__Basidiomycota;c__Malasseziomycetes;o__Malasseziales;f__Malasseziaceae;g__Malassezia",
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Cordycipitaceae",
                               "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Cantharellales",
                               "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Phaffomycetaceae;g__Starmera",
                               "k__Fungi;p__Basidiomycota;c__Malasseziomycetes;o__Malasseziales;f__Malasseziaceae;g__Malassezia",
                               "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Auriculariales;f__Exidiaceae",
                               "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus",
                               "k__Fungi;p__Ascomycota",
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Trichosphaeriales;f__Trichosphaeriaceae;g__Nigrospora",
                               "k__Fungi;p__Ascomycota",
                               "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Chaetothyriales;f__Cyphellophoraceae;g__Cyphellophora",
                               "k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Pleosporales;f__Didymellaceae;g__Epicoccum",
                               "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Saccharomycodaceae",
                               "k__Fungi;p__Ascomycota",
                               "k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Hymenochaetales;f__Tubulicrinaceae;g__Tubulicrinis",
                               "k__Fungi;p__Ascomycota;c__Eurotiomycetes;o__Eurotiales;f__Aspergillaceae;g__Aspergillus",
                               "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales",
                               "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Pichiaceae;g__Pichia",
                               "k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Sordariales;f__Chaetomiaceae;g__Mycothermus",
                               "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales",
                               "k__Fungi;p__Ascomycota",
                               "k__Fungi;p__Ascomycota;c__Saccharomycetes;o__Saccharomycetales;f__Trichomonascaceae")) %>% 
  separate_wider_delim(cols = Tax, delim = ";", too_few = "align_start",
                       names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), cols_remove = TRUE) %>% mutate(Species = NA) %>% 
  mutate(across(!OTU, ~str_remove(.x, "[:alpha:]__")))

# make new phylo
psnewname <- rarephy %>% tax_select(torem, deselect = TRUE)
newtt <- psnewname %>% tt_get() %>% as.data.frame() %>% rownames_to_column(var = "OTU") %>% filter(!OTU %in% torep) %>% 
  full_join(newnames) %>% column_to_rownames(var = "OTU") %>% as.matrix()
tax_table(psnewname) <- newtt

# save
save(psnewname, file = "data/phylo_ITS_resolvedNA.RData")

## ---- collapse ----

# get function
source("R/fx_myCollapse.R")

# collapse all ranks and remove redundant taxa (e.g.; a Genus with only one OTU in it is redundant with the OTU and should be removed)
allCol <- myCollapse(phylo = psnewname, tax_levels = c("OTU", "Genus", "Family", "Order", "Class", "Phylum"))

## collapse at only genus and above
#allColg <- myCollapse(phylo = rarephy,
#                 tax_levels = c("Genus", "Family", "Order", "Class", "Phylum"))


## ---- filter for prevalence ----

# remove any collapse taxa not present in at least 30% of individuals
pa <- decostand(allCol %>% 
                  column_to_rownames(var = "Sample"), method = "pa")  
csum <- colSums(pa)
tokeep <- names(which(csum > (nrow(allCol) * 0.30)))
allColfilt <- allCol %>% 
  select(Sample, all_of(tokeep)) %>% 
  column_to_rownames(var = "Sample")


## ---- normalize: CLR ----

# re-phyloseq for transformation
rephy <- phyloseq(
  otu_table(allColfilt, taxa_are_rows = FALSE),
  tax_table(data.frame(row.names = names(allColfilt), taxaname = names(allColfilt)) %>% as.matrix())
) %>% 
  tax_transform("clr", zero_replace = "halfmin")
# get data frame
allclr <- rephy %>% otu_get() %>% data.frame()


## ---- format for GWAS and sort ----

## the pheno file has to have the phenotype value in the 3rd column behind family and within-family IDs 
# can use the --pheno-name modifier to select a phenotype column by title or --mpheno to select numeric column, or use --all-pheno for all phenotypes in the pheno file
# in plink2, can remove FID and just have IID column


# wrangle; CLR
datCLR <- allclr %>% 
  rownames_to_column(var = "IID") %>% 
  mutate(FID = IID) %>% 
  relocate(FID, IID)


## this has to match the order of the .fam file
# sort and order
fam$IID <- as.character(fam$IID)
datCLR <- datCLR[order(match(datCLR$IID, fam$IID)), ]


# write to file
write.table(datCLR, "private/pheno_ITS.txt", sep = "\t", row.names = FALSE, quote = FALSE)
