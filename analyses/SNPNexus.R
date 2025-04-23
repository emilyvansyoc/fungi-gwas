## nov 2024 SNPNexus analysis
# EVS 12/2024

library(tidyverse)

### these results are from running SNPNexus in the web browser
# get results
df <- read.table("data/SNPNexus_results/allSNPNov2024_pervariant.tsv", sep = "\t", header = TRUE)

# get genomic coordinates with mine
coords <- read.table("data/SNPNexus_results/gen_coords_allSNPNov2024.txt", sep = "\t", header = TRUE)

# get my genomic coords
load("data/Nov2024_sigs_SNPandStructural.RData")

# match
both <- coords %>% 
  left_join(newsigs, by = c("Variation.ID" = "ID"))
match <- both %>% 
  filter(REF.Allele == REF & ALT.Allele..IUPAC. == ALT)

# get no matches
nom <- both %>% 
  anti_join(match) %>% 
  dplyr::select(Variation.ID, dbSNP, Chromosome, Position, REF.Allele, ALT.Allele..IUPAC., X.CHROM, POS, REF, ALT)

# make sure chromosomes match
nom %>% filter(!X.CHROM == Chromosome) # OK, all match (good!)

# get matching chromosomes and positions
nom1 <- nom %>% 
  filter(X.CHROM == Chromosome & Position == POS)

# get mismatched positions
nom2 <- nom %>% anti_join(nom1) #OK; thse are all structural variants that differ by 1 bp in position

# get matching reference alleles
nom %>% filter(REF.Allele == REF) %>% View()

## ---- get overlapped genes ----

# get data
genes <- read.table("data/SNPNexus_results/near_gens_allSNPNov2024.txt", sep = "\t", header = TRUE)

# get overlapped genes
unique(genes$Overlapped.Gene) # NAV2 is new 

# get overlapped genes
over <- genes %>% dplyr::filter(!Overlapped.Gene == "None") # 68

unique(over$Overlapped.Gene)
unique(over$Annotation) # none are exonic or otherwise coding
unique(over$Variation.ID) # 68 SNPs
unique(over$Chromosome) # four chromosomes; 1, 4, 11, 16

# double check that the genomic coordinates of these match with ours
coords1 <- coords %>% filter(Variation.ID %in% over$Variation.ID) %>% 
  left_join(newsigs, by = c("Variation.ID" = "ID"))
# nom
nom1 <- coords1 %>% filter(X.CHROM == Chromosome & Position == POS) # 62
nom2 <- coords1 %>% anti_join(nom1) # these are all structural
nom3 <- coords1 %>% filter(REF.Allele == REF & ALT.Allele..IUPAC. == ALT)
nom4 <- coords1 %>% anti_join(nom3) %>% relocate(Variation.ID, REF.Allele, REF, ALT.Allele..IUPAC., ALT)
nom5 <- coords1 %>% filter(REF.Allele == REF)

### investigate genes
genes %>% filter(!Overlapped.Gene == "None") %>% group_by(Variation.ID) %>% count() %>% arrange(desc(n))  #none are duplicated

## collapse these in the style of Rob's paper
genecol <- genes %>% 
  group_by(Variation.ID, Overlapped.Gene) %>% 
  summarize(anns = paste(Annotation, sep = ",")) %>% 
  filter(!Overlapped.Gene == "None")
length(unique(genecol$Variation.ID))
write.table(genecol, file = "data/results_SNPnexus_overlappedgenes_annotations.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# get fungal taxa associated with these genes
#load("R/gwas-output/Nov2024/Nov2024_sigs_SNPandStructural.RData")
genetax <- genecol %>% 
  left_join(newsigs, by = c("Variation.ID" = "ID"))
genetax %>% group_by(Overlapped.Gene, taxa) %>% count()

## ---- get ensembl ----

ens <- read.table("data/SNPNexus_results/ensembl_allSNPNov2024.txt", sep = "\t", header = TRUE)
unique(ens$Predicted.Function) # all are noncoding
