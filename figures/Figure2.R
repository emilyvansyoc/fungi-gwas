# Figure 2

library(tidyverse)
library(rstatix)
library(phyloseq)
library(microViz)
library(ggpubr)
library(vegan)
library(ggtext)
library(phyloseq)
library(microViz)
library(scales)

## ---- get fungal data ----

### get fungal data
load("data/phylo_ITS_resolvedNA.RData")
psf <- psnewname
# collapse and process via GWAS 
# get function
source("analyses/fx_myCollapse.R")

# collapse all ranks and remove redundant taxa (e.g.; a Genus with only one OTU in it is redundant with the OTU and should be removed)
allCol <- myCollapse(phylo = psnewname, tax_levels = c("OTU", "Genus", "Family", "Order", "Class", "Phylum"))


# remove any collapse taxa not present in at least 30% of individuals
pa <- decostand(allCol %>% 
                  column_to_rownames(var = "Sample"), method = "pa")  
csum <- colSums(pa)
tokeep <- names(which(csum > (nrow(allCol) * 0.30)))
allColfilt <- allCol %>% 
  select(Sample, all_of(tokeep)) %>% 
  column_to_rownames(var = "Sample")
# re-phyloseq for transformation
rephy <- phyloseq(
  otu_table(allColfilt, taxa_are_rows = FALSE),
  tax_table(data.frame(row.names = names(allColfilt), taxaname = names(allColfilt)) %>% as.matrix())
) %>% 
  tax_transform("clr", zero_replace = "halfmin")
# get data frame
allclr <- rephy %>% otu_get() %>% data.frame()

# get Kazachstania
kz <- allclr %>% dplyr::select(Genus_Kazachstania) %>% rownames_to_column(var = "IID") %>% mutate(IID = as.numeric(IID))


### ---- get gtes ----


load("data/Nov2024_GTEx_fullFAVResults_v10.RData")
### here, "slope" is equivalent to the GTEx portal "NES" effect size 
# so a positive slope is higher expression in the alt allele
# and a negative slope is higher expression in the ref allele

## get list of mycAVs to pull genotypes from PLINK
eqtls <- unique(full$rs_id_dbSNP151_GRCh38p7)
#write.table(eqtls, file = "data/Nov2024_GTEx_eQTLs.txt", row.names = FALSE, quote = FALSE)

# this is easier to do with separating single character from multi character alleles (n=78 of 82)
sing <- full %>% filter(nchar(ref) == 1 & nchar(alt) == 1)
eqtls.sing <- unique(sing$rs_id_dbSNP151_GRCh38p7)
#write.table(eqtls.sing, file = "data/Nov2024_GTEx_eQTLs_singlealleles.txt", row.names = FALSE, quote = FALSE)
mult <- full %>% anti_join(sing)
eqtls.mult <- unique(mult$rs_id_dbSNP151_GRCh38p7)
#write.table(eqtls.mult, file = "data/Nov2024_GTEx_eQTLs_multiplealleles.txt", row.names = FALSE, quote = FALSE)

## ---- v10 plots Jan 13 ----

# get genotypes
# headers: chromosome code, variant ID, position (can use dummy value of '0'), and base pair coordinates
map <- read.table("restricted/genos_Nov2024_GTEx_eQTLs.map", sep = "\t", header = FALSE, quote = "")
names(map) <- c("X.CHROM", "ID", "notneeded", "POS")
map1 <- map %>% dplyr::select(-notneeded) %>% 
  mutate(snpname = paste0("SNP", seq(1:nrow(map))))
# ped file headers; first 6 columns are FAM file, then alleles for first variant in .map file, etc
ped <- read.table("restricted/genos_Nov2024_GTEx_eQTLs.ped",  header = FALSE)
# remove extra columns
ped <- ped[,-c(2:6)]
map.ids <- unique(map$ID)

# read in with trio package
library(trio)
peddf <- read.pedfile("restricted/genos_Nov2024_GTEx_eQTLs.ped")
pedv <- peddf %>% 
  dplyr::select(-c(pid, fatid, motid, sex, affected)) %>% 
  dplyr::rename(IID = famid) %>% 
  pivot_longer(!IID, names_to = "snpID", values_to = "genotype") %>% 
  mutate(snpname = sapply(str_split(snpID, "\\."), `[`, 1)) %>% left_join(map1) %>% 
  dplyr::select(IID, ID, genotype)
# collapse
gens <- pedv %>% 
  group_by(IID, ID) %>% 
  mutate(geno.full = paste0(genotype, collapse = "")) %>% 
  dplyr::select(IID, ID, geno.full) %>% distinct()

### determine which are ref, alt, and heterozygous
genos <- gens %>% 
  left_join(full %>% dplyr::select(rs_id_dbSNP151_GRCh38p7, ref, alt) %>% distinct(), by = c("ID" = "rs_id_dbSNP151_GRCh38p7")) %>% 
  mutate(geno1 = sapply(substr(geno.full, 1, 1), `[`, 1),
         geno2 = sapply(substr(geno.full, 2, 2), `[`, 1)) %>% 
  mutate(geno.def = case_when(
    geno1 == ref & geno2 == ref ~ "REF",
    geno1 == alt & geno2 == alt ~ "ALT"
  )) %>% 
  mutate(geno.def = replace_na(geno.def, "Heterozy.")) %>% 
  mutate(geno.def = factor(geno.def, ordered = TRUE, levels = c("REF", "Heterozy.", "ALT")))

# add fungal rel abundance info
allclrv <- allclr %>% 
  rownames_to_column(var = "IID") %>% 
  pivot_longer(!IID, names_to = "taxa", values_to = "CLR")

# get eqtls with some wonky wrangling to fix issue with concatenating candida and saccharomycetales
eqtls <- full %>% 
  mutate(taxa1 = if_else(str_detect(taxa, "Candida"), "Genus_Candida", taxa)) %>% 
  rbind(full %>% filter(str_detect(taxa, "Candida")) %>% mutate(taxa1 = "Family_Saccharomycetales_fam_Incertae_sedis")) %>% 
  group_by(rs_id_dbSNP151_GRCh38p7, taxa1, gene_name, slope) %>% count() %>% dplyr::select(-n) %>% dplyr::rename(taxa = taxa1)

# add
plotdf <- genos %>% 
  left_join(eqtls, by = c("ID" = "rs_id_dbSNP151_GRCh38p7"), 
            relationship = "many-to-many") %>% 
  left_join(allclrv %>% filter(taxa %in% unique(eqtls$taxa)) %>% mutate(IID = as.numeric(IID)) %>% filter(IID %in% genos$IID)) %>% 
  # match taxa info (only 125 individuals)
  drop_na() %>% 
  mutate(posslope = if_else(sign(slope) == 1, "pos", "neg")) %>% 
  mutate(taxgene = paste(taxa, gene_name, sep = "-")) %>% 
  # make HHIP, HHIP-AS1, and OTUD4 into one "block"
  mutate(gene_cat = if_else(gene_name %in% "OTUD4", "HHIP, HHIP-AS1, OTUD4", gene_name)) %>% 
  filter(gene_cat %in% c("HHIP, HHIP-AS1, OTUD4", "CDH13", "GRIA1", "KRT8P33", "NAV2"))
plotdf <- plotdf %>% 
  mutate(taxa_genecat = paste(taxa, gene_cat, sep = "-"))

### sanity check: make sure the Candida genus and Saccharomycetales family are not completely identical
allclrv %>% filter(taxa %in% c("Genus_Candida", "Family_Saccharomycetales_fam_Incertae_sedis")) %>% pivot_wider(names_from = taxa, values_from = CLR) %>% mutate(match = if_else(Genus_Candida == Family_Saccharomycetales_fam_Incertae_sedis, "match", "nomatch")) %>% group_by(match) %>% count() ### only 15 people differ in the the abundance of these taxa; very similar 

###### DO STATS AND CONSTRUCT PLOTS #######

# need to know if fungal abundance varies with genotype? isn't that already determined in the GWAS?

ggplot(plotdf, aes(x = geno.def, y = CLR, fill = posslope)) +
  geom_boxplot() +
  #geom_point(position = position_jitter(width = 0.1, seed = 123)) +
  facet_wrap(~taxa_genecat, scales = "free") +
  theme_pubr()

### to circumvent the confusion over the directionality, make a category called "inc gene exp allele" instead of using ref/ale
### here, "slope" is equivalent to the GTEx portal "NES" effect size 
# so a positive slope is higher expression in the alt allele
# and a negative slope is higher expression in the ref allele

plotdf <- plotdf %>% 
  mutate(allele_cat = case_when(
    posslope == "pos" & geno.def == "ALT" ~ "Increased \nexpression \nalleles",
    posslope == "pos" & geno.def == "REF" ~ "Decreased \nexpression \nalleles",
    posslope == "neg" & geno.def == "ALT" ~ "Decreased \nexpression \nalleles",
    posslope == "neg" & geno.def == "REF" ~ "Increased \nexpression \nalleles",
    geno.def == "Heterozy." ~ "Heterozy."
  )) %>% 
  # pretty-fy category names
  mutate(pretty_cat = case_when(
    taxa_genecat %in% "Genus_Kazachstania-HHIP, HHIP-AS1, OTUD4" ~ "**B.** Kazachstania (*HHIP*, *HHIP-AS1*, *OTUD4*\\*)",
    str_detect(taxa_genecat, "Cand") ~ "**D.** Candida (*KRT8P33*, skin)",
    str_detect(taxa_genecat, "Sacc") ~ "**E.** Saccharomycetales (*KRT8P33*, skin)",
    str_detect(taxa_genecat, "Capno") ~ "**F.** Capnodiales (*GRIA1*, ganglia)",
    str_detect(taxa_genecat, "Asp") ~ "**C.** Aspergillaceae (*NAV2*, liver)",
    taxa_genecat %in% "Genus_Kazachstania-CDH13" ~ "**A.** Kazachstania (*CDH13*, arteries)"
  )) %>% 
  mutate(pretty_cat = factor(pretty_cat, ordered = TRUE, levels = c(
    "**A.** Kazachstania (*CDH13*, arteries)", "**B.** Kazachstania (*HHIP*, *HHIP-AS1*, *OTUD4*\\*)",
    "**C.** Aspergillaceae (*NAV2*, liver)", "**D.** Candida (*KRT8P33*, skin)", "**E.** Saccharomycetales (*KRT8P33*, skin)",
    "**F.** Capnodiales (*GRIA1*, ganglia)"
  )))

# run stats to add to the plot
stats <- plotdf %>%
  dplyr::group_by(pretty_cat) %>% 
  rstatix::dunn_test(CLR ~ allele_cat) %>% 
  add_significance() %>% 
  add_xy_position() %>% 
  mutate(ast = if_else(p.adj < 0.05, "*", "ns")) 

## ---- Feb 2025 add colors ----

# get the default ggplot2 colors from manhattan plots
cols <- hue_pal()(22)
names(cols) <- seq(1:22)

# add chromosome colors to points
col.names <- unique(plotdf1$pretty_cat)
new.names <- c(
  "A. <b style='color:#7C96FF'>Kazachstania</b> (*CDH13*, arteries)",
  "B. <b style='color:#C79800'>Kazachstania</b> (*HHIP*, *HHIP-AS1*, *OTUD4*\\*)",
  "C. <b style='color:#00C1A7'>Aspergillaceae</b> (*NAV2*, liver)",
  "D. <b style='color:#AEA200'>Candida</b> (*KRT8P33*, skin)",
  "E. <b style='color:#AEA200'>Saccharomycetales</b> (*KRT8P33*, skin)",
  "F. <b style='color:#AEA200'>Capnodiales</b> (*GRIA1*, ganglia)"
)
plotdf1 <- plotdf %>% 
  ungroup() %>% 
  mutate(newname = case_when(
    str_detect(pretty_cat, "A\\.") ~ new.names[1],
    str_detect(pretty_cat, "B\\.") ~ new.names[2],
    str_detect(pretty_cat, "C\\.") ~ new.names[3],
    str_detect(pretty_cat, "D\\.") ~ new.names[4],
    str_detect(pretty_cat, "E\\.") ~ new.names[5],
    str_detect(pretty_cat, "F\\.") ~ new.names[6]
  )) %>% 
  left_join(coords %>% dplyr::select(X.CHROM, ID, taxa)) %>% 
  left_join(data.frame(X.CHROM = as.numeric(names(cols)), hues = unname(cols)))

stats1 <- plotdf1 %>%
  dplyr::group_by(newname) %>% 
  rstatix::dunn_test(CLR ~ allele_cat) %>% 
  add_significance() %>% 
  add_xy_position() %>% 
  mutate(ast = if_else(p.adj < 0.05, "*", "ns")) 

# plot
ggplot(plotdf1, aes(x = allele_cat, y = CLR)) +
  geom_boxplot(size = 1, aes(fill = factor(X.CHROM))) +#, fill = "lightgrey") +
  #geom_point(aes(color = ID), position = position_jitter(width = 0.1, seed = 123), alpha = 0.2) +
  stat_pvalue_manual(stats1, label = "ast", hide.ns = TRUE, size = 6, bracket.size = 0.6, tip.length = 0) +
  labs(x = "**Tissue-specific gene expression alleles**", y = "**Fungal relative abundance**") +
  facet_wrap(~newname, scales = "free_x") +
  theme_pubr() +
  theme(strip.text = element_markdown(size = 12, hjust = 0),
        strip.background = element_rect(fill = "lightblue", color = NA),
        legend.position = "bottom",
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown()) +
  scale_fill_manual(values = cols, name = "Chromosome")

# save
ggsave(filename = "figures/eQTLs_fungi_color.png", dpi = 300, height = 8, width = 11, units = "in")
