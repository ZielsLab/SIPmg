library(tidyverse)
library(HTSSIP)
library(phyloseq)


# Pooling bin mapping statistics function
pool_bin_stat <- function(cov_tibble, bin_tibble, operator = "mean"){
  # DESCRIPTION: pool_bin_stat calculates per bin mean or median scaffold
  # value (coverage or mapped read numbers).
  #
  # INPUT: 
  #   cov_tibble - a tibble of scaffold values (numeric) that need to be pooled, 
  # column names should contain sample names, additional column "Feature" 
  # must be present and should contain scaffold IDs.
  #   bin_tibble - a tibble with two columns: first column with name "Feature" 
  # should contain scaffold IDs (same IDs used in cov_tibble), second column
  # with name "Bins" should contain bin IDs for each scaffold.
  #   operator - "mean" or "median" value, identifies pooling strategy.
  #
  # OUTPUT:
  #   a tibble with first column "Feature" that contains bin IDs, and the rest
  # of the columns represent samples with bins' pooled values.
  
  # setting the function to pool values based on "operator" parameter value
  dispatch_if <- function(operator, x){
    if(operator == "mean"){
      return(mean(x))
    } else if(operator == "median"){
      return(median(x))
    } else{
      stop("'operator' value must be either 'mean' or 'median' \n")
      break
    }
  }
  
  # join feature and bin tibbles based on "Feature" column and pooling values:  
  out_tibble <- inner_join(cov_tibble,bin_tibble, by="Feature") %>% 
    select(-Feature) %>% # removing scaffold names
    group_by(Bin) %>% # groupping based on assigned bins
    summarise_all(list(dispatch_if), operator = operator) # pooling values
  
  colnames(out_tibble)[1] <- "Feature" # removing scaffold names
  
  return(out_tibble)
}



