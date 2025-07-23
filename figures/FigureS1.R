### independent SNPs; results from all variants
# EVS 1/2024

library(tidyverse)
library(gridExtra)
library(ggpubr)
library(qqman)
library(phyloseq)
library(microViz)
library(cowplot)


# get effect sizes (too big to store in Github)
load("restricted/effectsizes_allSNP_allvariants_alltaxa.RData")



### ---- create and save QQ plots in a grid ----


mytax <- unique(alldat$taxa)
myplotlist <- c()

for(i in 1:length(mytax)) {
  
  # subset
  sub <- alldat %>% filter(taxa == mytax[i])
  
  # get observed and expected P values
  myp <- data.frame(
    observed = -log10(sort(sub$UNADJ)),
    expected = -log10(ppoints(nrow(sub)))
  )
  # get title
  mytitle <- str_remove_all(str_remove_all(mytax[i],
                                           "\\.glm\\.linear\\.adjusted"), "onlysnp\\.")
  mytitleFancy <- if_else(str_detect(mytitle, "OTU"),
                          # if OTU;
                          paste0(str_split(mytitle, "_")[[1]][1], ": ",
                                 str_split(mytitle, "_")[[1]][3]),
                          paste0(str_split(mytitle, "_")[[1]][1], ": ",
                         str_split(mytitle, "_")[[1]][4])
  )
  
  # plot with ggpubr
  myplot <- ggscatter(myp, x = "expected", y = "observed", alpha = 0.7,
            add = "reg.line", add.params = list(color = "red", linetype = "dashed", alpha = 0.8),
            xlab = "-Log10(Expected)", ylab = "-Log10(Observed)",
            title = mytitleFancy)
  
  # save to list
  myplotlist[[i]] <- myplot
  
  
}

### arrange and draw
gs <- ggarrange(plotlist = myplotlist, ncol = 6, nrow = 6)

# save
ggsave(plot = gs, dpi = 600, filename = "figures/allQQs.png",
       width = 20, height = 24, units = "in")
