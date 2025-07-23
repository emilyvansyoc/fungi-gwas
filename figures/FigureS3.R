# sunburst plot, updated 1/2025 with updated FAVs
# EVS 1/2025

library(plotly)
library(tidyverse)
library(ggpubr)

## ---- plot: sunburst of SNP Nexus gene annotations ----

# get data
# from Ensembl database
ens <- read.table("data/ensembl_allSNPNov2024.txt", sep = "\t", header = TRUE, comment.char = "")

# summarize
enssum <- ens %>% 
  # wrangle
  mutate(Symbol = if_else(Symbol %in% c("AC099506.1",
                                        "RN7SL134P"),
                          "CDH13", Symbol)) %>% 
  mutate(func = case_when(
    Predicted.Function %in% "intronic" ~ "Intronic",
    Predicted.Function %in% "non-coding intronic" ~ "Intronic",
    Predicted.Function %in% "non-coding" ~ "Intronic",
    Predicted.Function %in% "3downstream" ~ "Downstream 3'",
    Predicted.Function %in% "5upstream" ~ "Upstream 5'",
    Predicted.Function %in% "3utr" ~ "3' UTR"
  )) %>% 
  group_by(Symbol, func) %>% count()

tots <- enssum %>% ungroup() %>% group_by(Symbol) %>% summarize(tot = sum(n))
# build dataframe by hand
# add second hierarchy
mydf <- data.frame(
  ids = unique(enssum$Symbol),
  labels = unique(enssum$Symbol),
  parents = "",
  values = tots$tot
)
my2 <- data.frame(
  labels = enssum$func,
  values = enssum$n,
  parents = enssum$Symbol,
  ids = paste0(enssum$Symbol, " - ", enssum$func)
)
mydf <- rbind(mydf, my2) 

### wrangle and remove smaller sections
enssum1 <- enssum %>% ungroup() %>% 
  mutate(gene = if_else(Symbol %in% c("CDH13", "ANAPC10", "PTPRC", "NAV2"), Symbol, "Other"),
         func1 = if_else(func == "Intronic", "Intronic", "Other")) %>% 
  group_by(gene, func1) %>% summarize(tot = sum(n))
tots <- enssum1 %>% group_by(gene) %>% summarize(totgene = sum(tot))

mydf <- data.frame(
  ids = unique(enssum1$gene),
  labels = unique(enssum1$gene),
  parents = "",
  values = tots$totgene
)
my2 <- data.frame(
  labels = enssum1$func1,
  values = enssum1$tot,
  parents = enssum1$gene,
  ids = paste0(enssum1$gene, " - ", enssum1$func1)
)
mydf <- rbind(mydf, my2) 


## plot
sunb <- plot_ly() %>% 
  add_trace(data = mydf,
            labels = ~labels,
            parents = ~parents, 
            ids = ~ids, 
            values = ~values,
            type = "sunburst",
            branchvalues = 'total')


# save rough plot for now
reticulate::use_miniconda('r-reticulate')
plotly::save_image(sunb, file = "figures/FigureS3", scale = 3)
