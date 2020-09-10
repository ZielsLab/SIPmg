library(tidyverse)
library(HTSSIP)
library(phyloseq)


# Pooling bin mapping statistics function
pool_bin_stat <- function(cov_tibble, operator = "mean"){
  # DESCRIPTION: pool_bin_stat calculates per bin mean or median scaffold
  # value (coverage or mapped read numbers).
  #
  # INPUT: 
  #   cov_tibble - a tibble of coverage values (numeric) that need to be pooled, 
  # column names should contain Sequence Id, Bin Id, Coverage 
  # 
  #   
  #   operator - "mean" or "median" value, identifies pooling strategy.
  # 
  # OUTPUT:
  #   a tibble with first column "Feature" that contains bin IDs, and the rest
  # of the columns represent samples with bins' pooled values.
  
  # setting the function to pool values based on "operator" parameter value
f_tibble = bin_tibble
newnames <- f_tibble %>% select(grep(pattern = "Bam", x = names(.))) %>% #subset to BAM containing columns
  unique(.) %>% slice(1) %>% unlist(., use.names=FALSE) #get list of unique BAM file names
oldnames <- f_tibble %>% select(grep(pattern = "Coverage", x = names(f_tibble))) %>% names(.)
f_tibble <- f_tibble %>% rename_at(vars(oldnames), ~ newnames) #replace old names with new names
f_tibble <- f_tibble %>% select(1, 2, contains(newnames)) #pull out coverage columns
if (operator == "median") {
f_tibble = f_tibble %>%
  select(-`Sequence Id`) %>% 
  group_by(`Bin Id`) %>%
  summarise_all(median)
} else {
  f_tibble = f_tibble %>%
    select(-`Sequence Id`) %>% 
    group_by(`Bin Id`) %>%
    summarise_all(mean)
}
f_tibble = f_tibble %>%
  rename(Feature = `Bin Id`)
return(f_tibble)
}



