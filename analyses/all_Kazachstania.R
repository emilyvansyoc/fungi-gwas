### supplementary plot for all ordinations of kaz data
# EVS 1/2025

library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)
library(ggtext)
library(vegan)
library(EcolUtils)
set.seed(234)

### ---- build function to make plots ----

# set colors
kazcols <- c("darkgrey", "blue")
names(kazcols) <- c("Not detected", "Detected")

# build function
myOrdPlot <- function(ps.object) {
  
  ordplot <- ps.object %>% ps_filter(sample_sums(ps.object) > 1000) %>% # remove low seq depth to remove confounding factor
    ps_mutate(kazpres1 = if_else(kazpres == "Present", "Detected", "Not detected")) %>% 
    tax_transform("clr") %>% 
    ord_calc() %>% ord_plot(color = "kazpres1", size = 4, alpha = 0.6, auto_caption = NA
                            #plot_taxa = 1:8, 
                            #tax_vec_length = 1.3, tax_lab_length = 1.35, 
                            #tax_lab_style = list(size = 4, type = "text", max_angle = 90, fontface = "bold.italic", justify = "side") 
    ) +
    #stat_ellipse(aes(color = kazpres)) +
    #coord_fixed(ratio = 1, clip = "off", xlim = c(-3.8, 3.8), ylim = c(-3.5, 3)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title =  element_markdown(),
          text = element_text(size = 16),
          panel.background = element_rect(fill = NA, color = NA),
          plot.background = element_rect(fill = NA, color = NA)) +
    scale_color_manual(values = kazcols, name = "*Kazachstania*") +
    ggside::geom_xsideboxplot(aes(fill = kazpres1), orientation = "y", show.legend = FALSE) +
    ggside::geom_ysideboxplot(aes(fill = kazpres1), orientation = "x", show.legend = FALSE) +
    ggside::scale_xsidey_discrete(labels = NULL) +
    ggside::scale_ysidex_discrete(labels = NULL) +
    scale_fill_manual(values = kazcols) +
    ggside::theme_ggside_void()
}

## ---- build another function to rarefy ----

myRare <- function(ps.object, seqdepth) {
  
  # get otu table in tabular
  otab <- ps.object %>% otu_get() %>% as.data.frame()
  
  # rarefy 
  otab1 <- otab[rowSums(otab) > seqdepth, ]
  rtab <- EcolUtils::rrarefy.perm(otab1, sample = seqdepth, n = 100, round.out = TRUE)
  
  # re-phyloseq
  rarephy <- phyloseq(
    otu_table(rtab, taxa_are_rows = FALSE),
    tax_table(ps.object),
    sample_data(ps.object)
  )
  
  # remove empty taxa
  rarephy <- rarephy %>% tax_filter(min_prevalence = 1, min_total_abundance = 1)
  
  ### finish wrangling
  rarephy <- rarephy %>% tax_glom("Genus")
  taxa_names(rarephy) <- rarephy@tax_table[,'Genus']
  rarephy <- rarephy %>% 
    ps_otu2samdat("Kazachstania") %>% 
    ps_mutate(kazpres = if_else(Kazachstania > 0, "Present", "Absent"))
  
  # return rarefied object
  return(rarephy)
  
}

## ---- get all data and build plots ----

####### PRJNA541487 #######
load("publicITS_phyloseqs/raw_phyloseq_PRJNA541487.RData")
psf_541487 <- ps %>%
  tax_select("Fungi", "Kingdom") %>% 
  tax_fix()
sample_data(psf_541487) <- data.frame(row.names = sample_names(psf_541487)) %>% mutate(fakedata = "fakedata")

## get just healthy people
meta_541487 <- readxl::read_xlsx("publicITS_phyloseqs/Sra_PRJNA54187.xlsx") %>% 
  column_to_rownames(var = "Run")
sample_data(psf_541487) <- meta_541487
psf_541487g.healthy <- psf_541487 %>% 
  ps_filter(str_detect(`Sample Name`, "HC")) # 50 samples (paper says 52)

# this has very deep sequencing depth; rarefy at 30K
psf_541487g.healthy.rare <- myRare(psf_541487g.healthy, seqdepth = 30000) # 614 taxa

# what is percentage of kaz presence/absence?
psf_541487g.healthy.rare %>% samdat_tbl() %>% group_by(kazpres) %>% count() # 85% prevalence

# get stats
psf_541487g.healthy.rare %>% dist_calc("bray") %>% dist_permanova(variables = "kazpres") # significanct

# get ordination
plot_541487 <- myOrdPlot(psf_541487g.healthy.rare) +
  geom_text(label = "P = 0.005", color = "red", x = 1.5, y = 3.5,
            check_overlap = TRUE, fontface = "italic", size = 6)


##### PRJNA703732 ####
load("publicITS_phyloseqs/raw_phyloseq_PRJNA703732.RData")
psf_703732 <- ps %>% tax_fix() %>% tax_select("Fungi", "Kingdom")
## get just healthy controls
meta_703732 <- read.table("publicITS_phyloseqs/Sra_PRJNA703732.csv", sep = ",", header = TRUE) %>% column_to_rownames(var = "Run")
sample_data(psf_703732) <- meta_703732
psf_703732g.healthy <- psf_703732 %>% 
  ps_filter(str_detect(`sra_title..run.`, "control")) # 18 controls; matches paper

# decent sampling depth
psf_703732.rare <- myRare(psf_703732g.healthy, seqdepth = 7000)

# get percentage of kaz pres/abs
psf_703732.rare %>% samdat_tbl() %>% group_by(kazpres) %>% count()
# present in only 1 person

# plot
plot_703732 <- myOrdPlot(psf_703732.rare)

##### PRJNA698272 #####
load("publicITS_phyloseqs/raw_phyloseq_PRJNA698272.RData")
psf_698272 <- ps %>% tax_fix() %>% tax_select("Fungi", "Kingdom")
### get just healthy people
meta698272 <- read.table("publicITS_phyloseqs/Sra_PRJNA698272.csv", sep = ",", header = TRUE) %>% column_to_rownames(var = "Run")
sample_data(psf_698272) <- meta698272
psf_698272g.healthy <- psf_698272 %>% 
  ps_filter(str_detect(liver_disord, "healthy")) # 16 samples (paper states 16 healthy controls)

# very low sequencing depth
psf_698272.rare <- myRare(psf_698272g.healthy, seqdepth = 2000)

# get kaz pres/abs
psf_698272.rare %>% samdat_tbl() %>% group_by(kazpres) %>% count()
# present in 3 people

# stats
psf_698272.rare %>% dist_calc("bray") %>% dist_permanova(variables = "kazpres") # not significant

# plot
plot_698272 <- myOrdPlot(psf_698272.rare)

###### PRJEB11419 (American Gut Project)
load("publicITS_phyloseqs/raw_phyloseq_PRJEB11419.RData")
psf_11419 <- ps %>% tax_fix() %>% tax_select("Fungi", "Kingdom")
sample_data(psf_11419) <- data.frame(row.names = sample_names(psf_11419)) %>% mutate(fakedata = "fakedata")

# very low seq depth in many samples; rarefy to 20K
psf_11419.rare <- myRare(psf_11419, seqdepth = 5000) # lose 235 samples (over half) at 20K

# get kaz pres/abs
psf_11419.rare %>% samdat_tbl() %>% group_by(kazpres) %>% count() # present in only 2 people

# plot
plot_11419 <- myOrdPlot(psf_11419.rare)

#### ---- add together plots ----

mymargin <- theme(plot.margin = margin(t = 0.5, b = 0.5, unit = "cm"))

plot_541487 <- plot_541487 + mymargin
plot_698272 <- plot_698272 + mymargin
plot_703732 <- plot_703732 + mymargin
plot_11419 <- plot_11419 + mymargin

ggarrange(plot_541487, #plot_610042, 
          plot_698272, plot_703732, plot_11419, ncol = 2, nrow = 2,
          common.legend= TRUE, legend = "bottom",
          labels = c(
            "A. PRJNA541487 (Chinese)", 
            "B. PRJNA698272 (German)", "C. PRJNA703732 (Belgian)", "D. PRJEB11419 (American)"
          ),
          font.label = list(size = 16),
          hjust = 0, vjust = 0) +
  theme(plot.margin = margin(t=1, l = 0.2, r = 0.2, b = 0.2, unit = "cm"),
        plot.background = element_rect(fill = "white", color = NA))
# save
ggsave(filename = "figures/kazplots_alldata.png", dpi = 300, height = 12, width = 16, units = "in")
