## FIGURE 1
# EVS 12/2024

library(tidyverse)
library(ggpubr)
library(scales)
library(cowplot)
library(ggtext)
library(rcartocolor)
#library(BiocManager)
#BiocManager::install("ensembldb")
#BiocManager::install("EnsDb.Hsapiens.v75")
#install.packages("locuszoomr")
library(locuszoomr)
library(EnsDb.Hsapiens.v75)

set.seed(234)

### ----- panel A: Manhattan ----

# load FAVs
load("data/Nov2024_sigs_SNPandStructural.RData")
favs <- newsigs
# make identifier for taxa and SNPID
favs <- favs %>% 
  mutate(idsig = paste(taxa, ID, sep = "_"))

# get all data for manhattan plot (not shared; ginormous)
load("restricted/effectsizes_allSNP_allvariants_alltaxa.RData")

#### subset this artificially to make the plot smaller and more manageable
alleffects <- alleffects %>% 
  mutate(taxa = str_remove_all(taxa, "allsnp_allvariants\\."),
         taxa = str_remove_all(taxa, "\\.glm\\.linear")) %>% 
  mutate(idsig = paste(taxa, ID, sep = "_")) %>% 
  mutate(issig = if_else(idsig %in% favs$idsig, TRUE, FALSE))
# get the sigs
sigeffects <- alleffects %>% dplyr::filter(issig == TRUE)
# make smaller to reduce RAM
subdat <- alleffects %>% 
  #dplyr::filter(issig == FALSE) %>% 
  slice_sample(n=10000000)
#rm(alleffects) # need this in the locusZoomPlots script

# make even smaller for the plot - remove significant and filter p value to make plot prettier
subdat1 <- subdat %>% 
  dplyr::filter(issig == FALSE) %>% 
  dplyr::filter(-log10(P) > 1.7) %>% 
  slice_sample(n=10000)

# add back the significant points
plotdat <- subdat1 %>% 
  rbind(sigeffects)

# get Manhattan plot coordinates for chromosome names
coords <- plotdat %>% 
  group_by(X.CHROM) %>% 
  summarize(chr_len = max(POS)) %>% 
  ungroup() %>% 
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>% 
  dplyr::select(-chr_len) %>% 
  left_join(plotdat) %>% 
  arrange(X.CHROM, POS) %>% 
  mutate(BPcum = POS + tot) %>% 
  # add study-wise and exploratory P values for FAVs
  left_join(favs %>% dplyr::select(idsig, expsig, studysig)) %>% 
  # add thresholds for each significance 
  mutate(is.studysig = if_else(studysig < 0.2, TRUE, FALSE),
         is.gensig = if_else(P < 5e-8, TRUE, FALSE),
         is.expsig = if_else(expsig < 0.05, TRUE, FALSE))

# get positions of X axis titles
axisdf <- coords %>% group_by(X.CHROM) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)

# pretty-fy taxa names and get location of taxa text
sigtaxa <- coords %>% 
  dplyr::filter(issig == TRUE) %>% 
  mutate(taxa.pretty = case_when(
    taxa %in% "Order_Pleosporales" ~ "Pleosporales",
    taxa %in% "Genus_Kazachstania" ~ "Kazachstania",
    taxa %in% "Family_Saccharomycetaceae" ~ "Saccharomycetaceae",
    taxa %in% "Family_Saccharomycetales_fam_Incertae_sedis" ~ "Saccharomycetales I.S.",
    taxa %in% "Genus_Candida" ~ "Candida*",
    taxa %in% "Order_Capnodiales" ~ "Capnodiales",
    taxa %in% "Genus_Aspergillus" ~ "Aspergillus",
    taxa %in% "Family_Aspergillaceae" ~ "Aspergillaceae",
    taxa %in% "Order_Saccharomycetales" ~ "Saccharomycetales"
  )) %>% 
  group_by(taxa.pretty, X.CHROM) %>% 
  summarize(x.pos = max(BPcum),
            minP = min(P),
            minPlog = -log10(minP)) %>% 
  arrange(x.pos) 
# fix the "block" on chromosome 5 by using letters
sigtaxa.letters <- sigtaxa %>% 
  dplyr::filter(!taxa.pretty %in% c("Saccharomycetales I.S.")) %>% 
  cbind(plot.letters = LETTERS[1:10]) %>% 
  #cbind(plot.letters = c("A", "B", "C", "D \nE \nF", "G", "H", "I", "J", "K")) %>% 
  # nudge E by hand
  #mutate(x.pos = if_else(plot.letters == "E", x.pos + 1e7, x.pos))
  # nudge vertical instead of horizontal
  mutate(letter.y = if_else(plot.letters == "D", 9.2, 8.55),
         letter.y = if_else(plot.letters == "E", 9, letter.y))

## build plot
manplot <- ggplot(coords, aes(x = BPcum, y = -log10(P))) +
  # add approximate lines for signifance 
  # study-wise
  geom_hline(yintercept = -log10(1e-8), color = "blue") +
  # genome-wise
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  # exploratory
  geom_hline(yintercept = -log10(6.25e-7), color = "black") +
  
  # add points
  geom_point(aes(color = as.factor(X.CHROM)), alpha = 0.6, size = 3) +
  
  # add text for taxa names
  geom_segment(data = sigtaxa.letters, aes(x = x.pos, y = letter.y-0.2, yend = minPlog + 0.1), color = "grey") +
  geom_text(data = sigtaxa.letters, aes(x = x.pos, y = letter.y-0.2, label = taxa.pretty, color = as.factor(X.CHROM)), angle = 80,  hjust = "left", vjust = "center", size = 4.5) +
  
  # add labels for significance lines
  annotate("text", x = 2.95e9, y = -log10(1e-8) + 0.2, label = "study-wide", color = "blue", hjust = "right", fontface = "italic") +
  annotate("text", x = 2.95e9, y = -log10(5e-8) + 0.2, label = "genome-wide", color = "red", hjust = "right", fontface = "italic") +
  annotate("text", x = 2.95e9, y = -log10(6.25e-7) + 0.2, label = "exploratory", color = "black", hjust = "right", fontface = "italic") +
  
  # add text legend
  #geom_label(x = 3e9, y = 6, vjust = "top", hjust = "left", label = "A. Pleosporales \nB. Kazachstania \nC. Saccharomycetaceae \nD. Candida & Saccharomycetales I.S. \nE. Capnodiales \nF. Aspergillus \nG. Capnodiales \nH. Aspergillaceae \nI. Kazachstania \nJ. Saccharomycetales", size = 4, label.r = unit(0.1, "lines"), label.size = unit(0, "mm")) +
  
  ## back matter
  #ylim(1.7, 10.7) +
  #scale_y_continuous(limits = c(1.7, 10.7), breaks = c(2, 4, 6, 8)) +
  scale_x_continuous(label = axisdf$X.CHROM, breaks = axisdf$center) +
  theme_pubr(base_size = 14) +
  labs(x = "Chromosome", y = "-log<sub>10</sub>(*P*)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        axis.title.y = element_markdown(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        plot.margin = margin(t=3, unit = "cm")) +
  coord_cartesian(clip = "off", ylim = c(1.7, 8))


## ---- panel B-E: ZOOM plots ----

# API token
mytoken <- "mytoken"

### make plotting function

myLocusPlot <- function(locus.object, gene.name, showLegend = TRUE) {
  
  # get linkage diseq stats
  locus.object <- link_LD(locus.object, token = mytoken)
  
  # build plot
  #p <- locus_ggplot(locus.object, highlight = gene.name, 
  # show only certain types of genes for simplicity
  #filter_gene_biotype = c("protein_coding", "pseudogene"),
  #border = FALSE, showExons = TRUE, pcutoff = 0,
  #index_snp = NULL, #cex.lab = 2, cex.text = 1, 
  #shape = "isFAV", shape_values = c(24, 21), 
  # legend_pos = if_else(showLegend == TRUE, "left", "NULL"), size = 2.5)
  
  ## build plot with ggscatter and gg_genenames instead of locus_ggplot
  p <- gg_scatter(locus.object, highlight = gene.name,
                  border = FALSE, showExons = TRUE, pcutoff = 0,
                  index_snp = NULL, legend_pos = if_else(showLegend == TRUE, "bottom", "NULL"), size = 2.5,
                  shape = "isFAV", shape_values = c(24, 21), ylab = "-log<sub>10</sub>(*P*)") +
    theme_pubr(base_size = 14,
               legend = if_else(showLegend == TRUE, "bottom", "none")) +
    theme(#axis.line.x = element_blank(),
      #axis.ticks.x = element_blank(),
      axis.title.x = element_text(face = "bold"),
      #axis.text.x = element_blank(),
      plot.margin = margin(b=0, unit = "pt"),
      axis.title.y = element_markdown(face = "bold"))
  p1 <- gg_addgenes(p, locus.object, border = FALSE, highlight = gene.name, highlight_col = "black", gene_col = "grey",
                    exon_col = "grey", exon_border = "grey",
                    filter_gene_biotype = c("protein_coding"), cex.axis = 1.5, cex.lab = 1.5, cex.text = 1, xticks = FALSE) +
    theme_pubr(base_size = 14) +
    theme(plot.margin = margin(t=0, unit = "pt"))
  
  
  # return plot
  
  if(showLegend == TRUE) {return(p)} else{return(p1)}
  
  
}

### prep for plots

# GWAS data needs to be in dataframe format using 'locus'

## get full set of SNPs for most accurate plot for each taxa
# set pathway to output files
path <- "restricted/updateITS_results_allsnp_allvariant/"

## plot 1: chr1, pleosporales, gene "PTPRC"
pleo <- read.table(list.files(path, pattern = "Pleosporales", full.names = TRUE)[1], header = TRUE, sep = "\t", comment.char = "")
pleo <- pleo %>% 
  mutate(isFAV = factor(if_else(ID %in% favs$ID[favs$taxa %in% "Order_Pleosporales"], "FAV", "NS")))
# make locus object
ptprc.locus <- locus(gene = "PTPRC", data = pleo, ens_db = "EnsDb.Hsapiens.v75",
                     chrom = "X.CHROM", pos = "POS", p = "P", labs = "ID")
summary(ptprc.locus)

## plot2: chr 4, kazachstania, gene "ANAPC10"
kaz <- read.table(list.files(path, pattern = "Kazachstania", full.names = TRUE)[1], header = TRUE, sep = "\t", comment.char = "")
kaz <- kaz %>% 
  mutate(isFAV = factor(if_else(ID %in% favs$ID[favs$taxa %in% "Genus_Kazachstania"], "FAV", "NS")))
#make locus object
anap.locus <- locus(gene = "ANAPC10", data = kaz, ens_db = "EnsDb.Hsapiens.v75",
                    chrom = "X.CHROM", pos = "POS", p = "P", labs = "ID" )
summary(anap.locus)


## plot3: chr 11, aspergillaceae family, gene "NAV2"
asp <- read.table(list.files(path, pattern = "Aspergillaceae", full.names = TRUE)[1], header = TRUE, sep = "\t", comment.char = "")
asp <- asp %>% 
  mutate(isFAV = factor(if_else(ID %in% favs$ID[favs$taxa %in% "Family_Aspergillaceae"], "FAV", "NS")))
nav.locus <- locus(gene = "NAV2", data = asp, ens_db = "EnsDb.Hsapiens.v75",
                   chrom = "X.CHROM", pos = "POS", p = "P", labs = "ID" )
summary(nav.locus)


## plot4: chr 16, kaz, gene "CDH13"
cdh.locus <- locus(gene = "CDH13", data = kaz, ens_db = "EnsDb.Hsapiens.v75",
                   chrom = "X.CHROM", pos = "POS", p = "P", labs = "ID")
summary(cdh.locus)


### MAKE PLOTS

# make plot
chr1.plot <- myLocusPlot(locus.object = ptprc.locus, gene.name = "PTPRC", showLegend = FALSE)
# make plot
chr4.plot <- myLocusPlot(locus.object = anap.locus, gene.name = "ANAPC10", showLegend = FALSE)
# make plot
chr11.plot <- myLocusPlot(locus.object = nav.locus, gene.name = "NAV2", showLegend = FALSE)
# make plot
chr16.plot <- myLocusPlot(locus.object = cdh.locus, gene.name = "CDH13", showLegend = FALSE)

### ADD TOGETHER

allzooms <- ggarrange(chr1.plot, chr4.plot, chr11.plot, chr16.plot, ncol = 2, nrow = 2, 
                      labels = c("B. Pleosporales", "C. Kazachstania", "D. Aspergillaceae", "E. Kazachstania"),
                      font.label = list(size = 14, face = "bold", family = "Arial"),
                      hjust = -0.3, vjust = -0.5) +
  theme(plot.margin = margin(t = 0.7, r = 0, b = 0.2, l = 0, unit = "cm"))


leg <- myLocusPlot(locus.object = cdh.locus, gene.name = "CDH13", showLegend = TRUE)
#ggsave(filename = "figures/zoomplot_legend.png", dpi = 300)

leg1 <- ggpubr::get_legend(leg)

## add all and save
wleg <- ggarrange(allzooms, leg1, ncol = 1, heights = c(1, 0.1)) +
  theme(plot.margin = margin(t=0.5, unit = "cm"))

# add manhattan plot
ggarrange(manplot, wleg, ncol = 1, nrow = 2, heights = c(0.8, 1),
          labels = c("A. Fungi-associated variants", ""), 
          font.label = list(size = 14, face = "bold", family = "Arial"), hjust = -0.15, vjust = -0.2) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(t=1.5, unit = "cm"))

# save
ggsave(filename = "figures/wlegend_zoom_manhattan.png", height = 15, width = 10, units = "in")

