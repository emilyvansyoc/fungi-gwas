### FIGURE 3
# EVS 12/2024

library(tidyverse)
library(ggpubr)
library(phyloseq)
library(microViz)
library(vegan)
library(ggtext)
library(rstatix)

### --- panel A: kaz and heart disease alleles  ----

# get genotypes
gen <- read.table("restricted/extractsnps.ped", header = FALSE) %>% 
  dplyr::select(V1, V7, V8, V9, V10)
names(gen) <- c("IID", "rs12149890_a1", "rs12149890_a2", "rs12929586_a1", "rs12929586_a2")

gen <- gen %>% 
  mutate(rs12149890_geno = paste0(rs12149890_a1, rs12149890_a2),
         rs12929586_geno = paste0(rs12929586_a1, rs12929586_a2))

# get kazachstania abundance
load("data/phylo_ITS_resolvedNA.RData")
psf <- psnewname

# subset gen 
gensub <- gen %>% 
  dplyr::filter(IID %in% sample_names(psf))

# get Kazach CLR
kz <- psf %>% tax_glom("Genus") %>% 
  tax_transform("clr", zero_replace = "halfmin") %>% 
  tax_select("Kazachstania") %>% ps_melt() %>% dplyr::select(Sample, Abundance) %>% 
  mutate(IID = as.numeric(Sample))

# get Kazach prevalence
kzp <- psf %>% tax_glom("Genus") %>% 
  tax_transform("pa") %>% tax_select("Kazachstania") %>% ps_melt() %>% dplyr::select(Sample, Pres = Abundance)%>% 
  mutate(IID = as.numeric(Sample)) %>% 
  full_join(gensub) %>% 
  mutate(fullgeno = paste(rs12149890_geno, rs12929586_geno, sep = "/"))

# summarize
kzpsum <- kzp %>% group_by(Pres, fullgeno) %>% count()


### get data
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
  dplyr::select(Sample, all_of(tokeep)) %>% 
  column_to_rownames(var = "Sample")

# re-phyloseq for transformation
rephy <- phyloseq(
  otu_table(allColfilt, taxa_are_rows = FALSE),
  tax_table(data.frame(row.names = names(allColfilt), taxaname = names(allColfilt)) %>% as.matrix())
) %>% 
  tax_transform("clr", zero_replace = "halfmin")
# get data frame
allclr <- rephy %>% otu_get() %>% data.frame()

# join and plot
forplot <- gensub %>% 
  #full_join(kz) %>% 
  full_join(allclr %>% dplyr::select(Genus_Kazachstania) %>% rownames_to_column(var = "IID") %>% mutate(IID = as.numeric(IID))) %>% 
  mutate(fullgeno = paste(rs12149890_geno, rs12929586_geno, sep = "/")) %>% 
  mutate(fullgeno = if_else(fullgeno == "TG/AC", "GT/CA", fullgeno)) %>% 
  mutate(fullgeno = factor(fullgeno, ordered = TRUE, levels = c("GG/CC", "GT/CA", "TT/AA")))

# get sigs
stats <- forplot %>% dunn_test(Genus_Kazachstania ~ fullgeno) %>% 
  add_significance() %>% 
  add_xy_position() %>% 
  mutate(y.position = y.position - 2)

### get fold change between homozygous and heterozygous
het <- summary(forplot$Genus_Kazachstania[forplot$fullgeno == "GT/CA"])
hom <- summary(forplot$Genus_Kazachstania[forplot$fullgeno == "GG/CC"])
hom.alt <- summary(forplot$Genus_Kazachstania[forplot$fullgeno == "TT/AA"])
hom/het # on average 2.6x lower in risk allele compared to heterozygous
hom/hom.alt #not useful due to small sample size


### color the x axis text for variant alleles
plotcol <- forplot %>% 
  mutate(a1col = case_when(
    fullgeno == "GG/CC" ~ "#FF0000",
    fullgeno == "GT/CA" ~ "#FF0000",
    fullgeno == "TT/AA" ~ "#000000"
  ),
  a2col = case_when(
    fullgeno == "GG/CC" ~ "#FF0000",
    fullgeno == "GT/CA" ~ "#000000",
    fullgeno == "TT/AA" ~ "#000000"
  ),
  rs12149890_altal = case_when(
    rs12149890_geno == "GG" ~ "G",
    rs12149890_geno == "TG" ~ "G",
    rs12149890_geno == "TT" ~ "T"
  ),
  rs12929586_altal = case_when(
    rs12929586_geno == "CC" ~ "C",
    rs12929586_geno == "AC" ~ "C",
    rs12929586_geno == "AA" ~ "A"
  ),
  rs12149890_refal = case_when(
    rs12149890_geno == "GG" ~ "G",
    rs12149890_geno == "TG" ~ "T",
    rs12149890_geno == "TT" ~ "T"
  ),
  rs12929586_refal = case_when(
    rs12929586_geno == "CC" ~ "C",
    rs12929586_geno == "AC" ~ "A",
    rs12929586_geno == "AA" ~ "A"
  )) %>% 
  mutate(
    collab = paste0("<b style='color:", a1col, "'>", rs12149890_altal, "</b><b style='color:", a2col, "'>", rs12149890_refal, "</b><b style='color:#000000'>/</b><b style='color:", a1col, "'>", rs12929586_altal, "</b><b style='color:", a2col, "'>", rs12929586_refal, "</b>")
  ) 

# order correctly
unique(plotcol$collab)
plotcol <- plotcol %>% mutate(collab = factor(collab, ordered = TRUE, levels = c("<b style='color:#FF0000'>G</b><b style='color:#FF0000'>G</b><b style='color:#000000'>/</b><b style='color:#FF0000'>C</b><b style='color:#FF0000'>C</b>", 
                                                                                 "<b style='color:#FF0000'>G</b><b style='color:#000000'>T</b><b style='color:#000000'>/</b><b style='color:#FF0000'>C</b><b style='color:#000000'>A</b>",
                                                                                 "<b style='color:#000000'>T</b><b style='color:#000000'>T</b><b style='color:#000000'>/</b><b style='color:#000000'>A</b><b style='color:#000000'>A</b>")))


## build plot
kazboxplot <- ggplot(plotcol, aes(x = collab, y = Genus_Kazachstania)) +
  geom_boxplot(outliers = FALSE, size = 1) +
  geom_point(color = "blue", alpha = 0.4, size = 4, position = position_jitter(width = 0.2, seed = 123)) +
  labs(x = "rs12149890 / rs12929586 genotypes", y = "*Kazachstania* abundance") +
  stat_pvalue_manual(stats, label = "p.adj.signif", hide.ns = TRUE, size = 8, bracket.size = 0.8, bracket.nudge.y = 0) +
  theme_pubr() +
  theme(axis.title.y = element_markdown(size = 14, face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.x = element_markdown(size = 16),
        text = element_text(size = 16))+
  coord_fixed(ratio = 0.15)


### ---- panel B: kaz prevalence and abundance ----

psg <- psf %>% tax_fix() %>% tax_glom("Genus")
taxa_names(psg) <- psg@tax_table[,"Genus"]

# get only the top 25 taxa
topn <- tax_top(psg, n = 15, rank = "Genus", by = "prev")

# get prevalence of all genera
kzp <- psg %>% 
  tax_select(topn, "Genus", strict_matches = TRUE) %>% 
  tax_transform("pa") %>% ps_melt() %>%
  dplyr::select(Sample, Pres = Abundance, OTU) %>% 
  # calculate total prevalence
  group_by(OTU) %>% summarize(tot = sum(Pres)) %>% mutate(prev = (tot/125)*100) %>% 
  ungroup() %>% 
  mutate(OTU = factor(OTU, ordered = TRUE, levels = rev(topn))) %>% 
  mutate(iskaz = if_else(OTU == "Kazachstania", "kaz", "not")) 

# get relative abundance of all genera
ra <- psg %>% 
  tax_select(topn, "Genus", strict_matches = TRUE) %>% 
  tax_transform("compositional") %>% 
  ps_melt() %>% dplyr::select(Sample, Abundance, OTU) %>% 
  mutate(OTU = factor(OTU, ordered = TRUE, levels = rev(topn))) %>% 
  mutate(iskaz = if_else(OTU == "Kazachstania", "kaz", "not")) 

## plot
# boxplot for RA
## super hacky way to change Kazachsatnia color to red
labs <- paste("<span style = 'color: ",
              ifelse(topn == "Kazachstania", "blue", "black"),
              ";'>",
              topn,
              "</span>", sep = "")
ra <- ra %>% 
  mutate(mylab = paste("<span style = 'color: ",
                       ifelse(iskaz == "kaz", "blue", "black"),
                       ";'>",
                       OTU,
                       "</span>", sep = "")) %>% 
  mutate(mylab = factor(mylab, ordered = TRUE, levels = rev(labs)))

## 1/9/25; swap order so prevalence is left and RA is right
pa <- ggplot(data = ra, aes(x = Abundance, y = mylab, color = iskaz)) +
  #geom_boxplot() + 
  geom_point(position = position_jitter(height = 0.1, seed = 123), alpha = 0.6) +
  theme_pubr() +
  xlab("Relative Abundance") + ylab("") +
  theme(legend.position = "none",
        axis.text.y = element_blank(),
        axis.title.x = element_text(face = "bold"),
        axis.ticks.y = element_blank(),
        text = element_text(size = 16),
        panel.background = element_rect(fill = NA, color = NA),
        plot.background = element_rect(fill = NA, color = NA)) +
  scale_color_manual(values = c("blue", "black"))

# barplot for prevalence
pb <- ggplot(data = kzp, aes(x = prev, y = OTU, fill = iskaz)) +
  geom_bar(stat = "identity") +
  #scale_x_reverse(limits = c(100, 0)) +
  xlim(c(0, 100)) +
  theme_pubr() +
  xlab("Prevalence (%)") + ylab("") +
  scale_fill_manual(values = c("blue", "lightgrey")) +
  #guides(y = "none", y.sec = "axis") + # move axis line to the right side
  
  theme(axis.text.y = element_markdown(face = "italic", size = 14),
        axis.title.x = element_text(face = "bold"),
        legend.position = "none",
        text = element_text(size = 16),
        panel.background = element_rect(fill = NA, color = NA),
        plot.background = element_rect(fill = NA, color = NA))

## arrange
prevfig <- ggarrange(pb, pa, ncol = 2, widths = c(1, 0.6)) +
  theme(panel.background = element_rect(fill = NA, color = NA))

### ---- panel C: ordination of kaz presence/absence ----

psg <- psg %>% ps_otu2samdat("Kazachstania") %>% 
  ps_mutate(kazpres = if_else(Kazachstania > 0, "Detected", "Not detected")) %>% 
  ps_mutate(seqdepth = sample_sums(psg))

# get stats
psg %>% dist_calc("bray") %>% dist_permanova(variables = c("kazpres", "seqdepth")) # seqdepth barely not sig
psg %>% ps_filter(seqdepth > 1300) %>% dist_calc("bray") %>% dist_permanova(variables = c("kazpres", "seqdepth"))

kazcols <- c("darkgrey", "blue")
names(kazcols) <- c("Not detected", "Detected")

## build ordination (use PCA to show taxa loadings)
ordplot <- psg %>%# ps_filter(sample_sums(psg) > 1000) %>% # remove low seq depth to remove confounding factor
  tax_transform("clr") %>% 
  ord_calc() %>% ord_plot(color = "kazpres", size = 4, alpha = 0.6, auto_caption = NA,
                          plot_taxa = 1:8, 
                          tax_vec_length = 0.95, tax_lab_length = 1, 
                          tax_lab_style = list(size = 4, type = "text", max_angle = 90, fontface = "bold.italic", 
                                               justify = "side") ) +
  #stat_ellipse(aes(color = kazpres)) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.2, 3), ylim = c(-3, 3)) +
  # add permanova stat in corner
  #ggtext::geom_richtext(x = 2.5, y = 2, label = "PERMANOVA <br>R<sup>2</sup> = 0.15, *P* = 0.001", size = 4, label.size = 0) +
  theme_pubr() +
  theme(legend.position = "bottom",
        legend.title =  element_markdown(),
        text = element_text(size = 16),
        panel.background = element_rect(fill = NA, color = NA),
        plot.background = element_rect(fill = NA, color = NA),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold")) +
  scale_color_manual(values = kazcols, name = "*Kazachstania*") +
  ggside::geom_xsideboxplot(aes(fill = kazpres), orientation = "y", show.legend = FALSE) +
  ggside::geom_ysideboxplot(aes(fill = kazpres), orientation = "x", show.legend = FALSE) +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  scale_fill_manual(values = kazcols) +
  ggside::theme_ggside_void()

##---- arrange and save ----



left <- ggarrange(kazboxplot, prevfig,  ncol = 1, heights = c(0.5, 1), labels = c("A.", "B."), hjust = -3,
                  font.label = list(size = 16, face = "bold")) 
ggarrange(left, ordplot, ncol = 2,  labels = c("", "C."), font.label = list(size = 16, face = "bold")) +
  theme(plot.background = element_rect(fill = "white", color = "white"),
        plot.margin = margin(t=0.5, unit = "cm"))
ggsave(filename = "figures/Nov2024_kazplot.png", dpi = 600)
