library(tidyverse)
library(ggpubr)

# Sequin global scaling function

scale_features <- function(f_tibble, sequin_meta, seq_dilution, log_trans){
  
  # DESCRIPTION: scale_features calculates global scaling factors for features (contigs or bins),
  # based on linear regression of sequin coverage. Options include log-transformations of 
  # coverage, as well as filtering features based on limit of detection. This function must be called
  # first, before the feature abundance table, feature detection table, and plots are retrieved.
  # 
  #  INPUT: 
  #   f_tibble -    can be either: 
  #       (1) a tibble with first column "Feature" that contains bin IDs, and the rest
  #       of the columns represent samples with bins' pooled values. Every sequin is also listed
  #       as a feature.
  #       (2) a tibble with first column "Feature" that contains contig IDs,and the rest
  #       of the columns represent samples with contig coverage values. Every sequin is also listed
  #       as a feature.
  #   sequin_meta -     tibble containing sequin names ("Feature column") and concentrations in 
  #       attamoles/uL ("Concentration") column.
  #   log_trans -       Boolean (TRUE or FALSE), should coverages and sequin concentrations 
  #       be log-scaled?
  #   seq_dilution -    tibble with first column "Sample" with *same sample names as in f_tibble*, and a
  #       second column "Dilution" showing ratio of sequins added to final sample volume 
  #       (e.g. a value of 0.01 for a dilution of 1 volume sequin to 99 volumes sample)
  #
  # OUTPUT:
  #   mag_tab -       a tibble with first column "Feature" that contains bin (or contig IDs), 
  #       and the rest of the columns represent samples with features' scaled abundances (attamoles/uL)
  #   mag_det -       a tibble with first column "Feature" that contains bin (or contig IDs), 
  #       and the rest of the columns represent samples with scaled abundance for those above LOD, or 
  #       NA if it is below detection.
  #   plots   -       a tibble of spike-in calibration curves for all samples
  #   scale_fac -     a master tibble with all of the intermediate values in above calculations
  

  # Retrieve sample names from feature tibble
    scale_fac <- tibble(Sample = names(f_tibble)) %>% filter(!grepl("Feature", Sample))
  
  # Merge dilution factors for samples, add log-scaling option
    scale_fac <- scale_fac %>% inner_join(seq_dilution, by = "Sample") %>% mutate(log_scale = log_trans)
    
  # Make coverage table for features
    scale_fac <- scale_fac %>% 
      mutate(
        cov_tab = map(Sample, ~ select(f_tibble, Feature, all_of(.))), # Make list of coverage tables for samples
        cov_tab = map(cov_tab, ~setNames(., c("Feature", "Coverage")))  #get rid of sample name in header
        ) %>%

      # merge sequins with coverage tables, remove them from mag/feature table
      mutate( 
        seq_cov = map(cov_tab, ~ inner_join(., sequins, by = "Feature")), 
        mag_cov = map(cov_tab, ~ anti_join(., sequins, by = "Feature"))
      )  %>%
      
      # scale sequin concentrations based on dilution factors
      mutate(
        seq_cov = map2(seq_cov, Dilution, ~ mutate(.x, Concentration = Concentration / .y))
      ) %>%
    
      # Determine groups of spike-in concentrations 
      mutate(
      seq_group = map(seq_cov, ~.x %>% group_by(Concentration) %>% tally(name="standards"))
      ) %>%
      
      # determine limit of detection of spike-ins, based on presence of 5 sequins per conc. 
      mutate(
        #determine lowest concentration where at least 1 sequins is detected 
        lod = map_dbl(seq_cov, ~ filter(., Coverage > 0) %>% group_by(., Concentration) %>% 
                        tally(name="detected") %>% summarise(Min = min(Concentration)) %>% pull(Min)),
        
        #create tibble comparing number of observed and theoretical spike ins
        seq_det = map2(seq_cov, seq_group, ~ filter(., Coverage > 0) %>% group_by(., Concentration) %>% 
                        tally(name="detected") %>% inner_join(.y, by = "Concentration")), 
        
        # determine difference between standards and observed spike ins
        seq_det = map(seq_det, ~ mutate(., diff = standards - detected)),
        seq_warning = map_int(seq_det, ~summarise(., Sum = sum(diff)) %>% pull(Sum)) #positive values give warning later
      ) %>%
      
      # perform linear regression on coverage vs conc., extract lm params, make plots
      mutate(
        seq_cov_filt = map(seq_cov, ~ filter(., Coverage > 0)), #remove zero coverage values before lm
                            
        fit = ifelse(log_scale == "TRUE" , # check log_trans input
                      map(seq_cov_filt, ~ lm(log10(Concentration) ~ log10(Coverage) , data = .)), #log lm if true
                      map(seq_cov_filt, ~ lm((Concentration) ~ (Coverage) , data = .)) # lm if false
                      ), 
        slope = map_dbl(fit, ~summary(.)$coef[2]), # get slope
        intercept = map_dbl(fit, ~summary(.)$coef[1]) # get intercept
        ) %>%
      
      #plot linear regressions
      mutate(
        plots = ifelse(log_scale == "TRUE" , # check log_trans input
                       map(seq_cov_filt, # log-scaled plot if true
                         ~ ggplot(data=. , aes(x=log10(Coverage), y= log10(Concentration))) + 
                           geom_point() + 
                           geom_smooth(method = "lm") + 
                           stat_regline_equation(label.x= -0.1, label.y = 3) + 
                           stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -0.1, label.y = 3.5) + 
                           xlab("Coverage (log[read depth])") + 
                           ylab("DNA Concentration (log[attamoles/uL])") + 
                           theme_bw() 
                            ),
                         map(seq_cov_filt, # non-scaled plot if true
                             ~ ggplot(data=. , aes(x=Coverage, y= Concentration)) + 
                               geom_point() + 
                               geom_smooth(method = "lm") + 
                               stat_regline_equation(label.x= 0, label.y = 1000) + 
                               stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0, label.y = 1500) + 
                               xlab("Coverage (read depth)") + 
                               ylab("DNA Concentration (attamoles/uL)") + 
                               theme_bw()
                             )
                         )
        ) %>%
      
      #flag MAGs below LOD, and scale MAGs by slope and intercept 
      mutate(
        # Scale MAGs based on linear regression
        mag_ab = ifelse(log_scale == "TRUE" , # check log_trans input
                        map2(mag_cov, slope, ~ mutate(.x, Concentration = log10(Coverage) * .y)), # y = mx (in log scale) if true
                        map2(mag_cov, slope, ~ mutate(.x, Concentration = Coverage * .y)) # y = mx if false
        ),
        mag_ab = map2(mag_ab, intercept, ~ mutate(.x, Concentration = Concentration + .y)), # +b
        mag_ab = ifelse(log_scale == "TRUE" , # check log_trans input
                        map(mag_ab, ~ mutate(.x, Concentration = 10^Concentration)), #convert back from log10 if true
                        mag_ab), # no change if false
        mag_ab = map(mag_ab, ~ select(., Feature, Concentration)), # drop Coverage column
        mag_ab = map2(mag_ab, Sample, ~ setnames(.x, 'Concentration', .y)), #put sample name in MAG table
        
        # Remove MAGs below LOD
        mag_det = mag_ab,
        mag_det = map2(mag_det, Sample, ~ setnames(.x, .y, 'Concentration')), #get sample name out of header for filter
        mag_det = map2(mag_det, lod, ~ filter(.x, Concentration > .y)),
        mag_det = map2(mag_det, Sample, ~ setnames(.x, 'Concentration', .y)) #change header back to sample
        )
    
    # compile feature abundance across samples
    mag_tab <- scale_fac$mag_ab %>% reduce(left_join, by="Feature")
    mag_names <- mag_tab[[1]]
    mag_tab <- data.matrix(mag_tab[,-1])
    rownames(mag_tab) <- mag_names
    
    # extract plots
    plots <- scale_fac %>% select(Sample, plots)
    
    # extract feature detection 
    mag_det <- scale_fac$mag_det %>% reduce(left_join, by="Feature")
    mag_dnames <- mag_det[[1]]
    mag_det <- data.matrix(mag_det[,-1])
    rownames(mag_det) <- mag_dnames
    
    # make list of results
    results <- list("mag_tab" = mag_tab, 
                    "mag_det" = mag_det, 
                    "plots" = plots, 
                    "scale_fac" = scale_fac)
    
    # return results
    return(results)

}

### Example usage #### 

f_tibble <- read_csv(file="mock_import/pool_bin_stat_w_seqins.csv")
sequins <- read_csv(file="mock_import/sequins_metadata.csv")
seq_dil <- tibble( Sample = c("S1", "S2", "S3", "S4"), Dilution = c(0.01, 0.01, 0.05, 0.05)) 
log_scale <- "TRUE"

mag_tab_scaled <- scale_features(f_tibble, sequins, seq_dil, log_scale)

