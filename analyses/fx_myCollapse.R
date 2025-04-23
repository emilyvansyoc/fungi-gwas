### function to collapse phyloseq object by taxonomic levels
# EVS 1/2024

library(microViz)
library(phyloseq)
library(tidyverse) 

# build function

myCollapse <- function(phylo, # input is a phyloseq object
                       tax_levels = c("OTU", "Genus", "Family", "Order", "Class", "Phylum")) # and tax levels to collapse at  
  { 
  
  # dataframe-ize
  df <- phylo %>%
    #tax_fix() %>% 
    ps_melt() %>% 
    dplyr::select(Sample, Abundance, all_of(tax_levels))
  

  ## collapse at each tax level
  outdf <- data.frame(Sample = unique(df$Sample))
  mylevels <- unique(tax_levels)
  
  for(i in 1:length(mylevels)) {
    
    # get redundant taxa levels that can be ignored
    if(i > 1) {
      torem <- df %>% dplyr::select(mylevels[i-1], mylevels[i]) %>% distinct() 
      names(torem) <- c("mysublevel", "myLevel")
      torem1 <- torem %>% drop_na(myLevel) %>% group_by(myLevel) %>% count() %>% 
        dplyr::filter(n == 1) %>% dplyr::select(myLevel) %>% ungroup()
    }
    
    # collapse
    col <- df %>% 
      dplyr::select(Sample, Abundance, mylevels[i]) 
    names(col) <- c("Sample", "Abundance", "myLevel")
    col <- col %>%  
      drop_na(myLevel) %>% 
      group_by(Sample, myLevel) %>% 
      dplyr::summarize(tot_tax = sum(Abundance)) 
    
    # remove redundant taxa
    if(i > 1) {
      colsub <- col %>% anti_join(torem1) %>% 
        # pivot wider
        pivot_wider(names_from = myLevel, values_from = tot_tax, names_prefix = paste0(mylevels[i], "_"))
    } else {
      colsub <- col %>% 
        # pivot wider
        pivot_wider(names_from = myLevel, values_from = tot_tax, names_prefix = paste0(mylevels[i], "_"))
    }
 
    # save
    outdf <- full_join(outdf, colsub)
    
  }
  
  
  ## replace any spaces with a dash
  names(outdf) <- gsub(" ", "-", names(outdf))
  
  
  # return joined dataset
  return(outdf)
  
}
