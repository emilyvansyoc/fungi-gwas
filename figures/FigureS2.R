### SNP linkage disequilibrium 
#EVS 5/1/2024, updated 1/15/2025

library(tidyverse)
library(ggpubr)
library(ggtext)


myids <- read.table("data/unique_FAVs.txt", sep = "\t", header = TRUE) 
# write out to get LD report in PLINK
#write.table(myids$ID, file = "data/unique_FAVs_plink.txt", sep = "\t", row.names = FALSE, quote = FALSE)

## ---- ld report ----
# get LD report
ld <- read.table("data/ld_report_snps_Nov2024.ld", header = TRUE, quote = "", comment.char = "")

# remove duplicates (fixed in plink2 but hasn't been implemented)
ld1 <- ld %>% filter(!SNP_A == SNP_B)

# get only the variants in our mycAVs (not sure why extras are added?)
ld2 <- ld1 %>% 
  filter(SNP_A %in% myids$ID & SNP_B %in% myids$ID)

# make horizontal
ldh <- ld2 %>% 
  dplyr::select(SNP_A, SNP_B, R2) %>% 
  pivot_wider(names_from = SNP_B, values_from = R2)

## ---- get ld for each SNP ----
# get all significant FAVs with their taxa
load("data/Nov2024_sigs_SNPandStructural.RData")

# this will only work for taxa associated with more than one FAV
newsigs %>% group_by(taxa) %>% count() %>% filter(n==1)

# add taxa info and make vertical
both <- newsigs %>% dplyr::select(taxa, ID, X.CHROM) %>% 
  left_join(ldh, by = c("ID" = "SNP_A"))%>% 
  pivot_longer(cols = starts_with("rs"), names_to = "SNP_B", values_to = "R2") %>% 
  drop_na(R2)

###  plot
all <- both %>% 
  mutate(Taxa = str_replace(taxa, "_", " ")) %>% 
  #mutate(Taxa = if_else(str_detect(Taxa, "Kaz"), "Genus *Kazachstania*", Taxa)) %>% 
  #mutate(Taxa = factor(Taxa, ordered = TRUE, levels = c(
  # "Order Pleosporales", "Order Saccharomycetales", "Order Capnodiales", "Genus *Kazachstania*"
  #))) %>% 
  mutate(Taxa = fct_reorder(Taxa, X.CHROM)) %>% 
  mutate(X.CHROM = factor(X.CHROM))

# make plot
ggdensity(all, x = "R2", fill = "X.CHROM", facet.by = "Taxa", scales = "free_y",
          #bins = 10, 
          xlab = "R<sup>2</sup>", ylab = "Number of FAVs",
          position = "stack") +
  guides(fill = guide_legend(title = "Chromosome")) +
  theme_pubr(base_size = 14, legend = "right") +
  theme(strip.text = element_markdown(size = 18),
        axis.title.x = element_markdown(),
        strip.background = element_rect(fill = "white"))
ggsave(filename = "figures/FigureS2.png", dpi = 300)
