library(tidyverse)
library(ggpubr)

# Sequin global scaling function

scale_features_ps <- function(f_tibble, sequin_meta, seq_dilution, log_trans, coe_of_variation){
  
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
  #       (2) a tibble as outputted by the program `checkm coverage`. If this is the
  #       input format, the optional function, pooling_functions.R must be run.
  #       `pooling_functions.R` parses the checkM coverage output to provide a tibble as described in option 1
  #       Please check `pooling_functions.R` for further details. Please check checkm documentation on the usage 
  #       for `checkm coverage` program
  #   sequin_meta -     tibble containing sequin names ("Feature column") and concentrations in 
  #       attamoles/uL ("Concentration") column.
  #   log_trans -       Boolean (TRUE or FALSE), should coverages and sequin concentrations 
  #       be log-scaled?
  #   seq_dilution -    tibble with first column "Sample" with *same sample names as in f_tibble*, and a
  #       second column "Dilution" showing ratio of sequins added to final sample volume 
  #       (e.g. a value of 0.01 for a dilution of 1 volume sequin to 99 volumes sample)
  #   threshold - Acceptable coefficient of variation for coverage and detection 
  #       (eg. 20 - for 20 % threshold of coefficient of variation)
  #       (Coverages above the threshold value will be flagged in the plots)
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
  scale_fac <- tibble(Sample = names(f_tibble) %>% str_subset(pattern = "Feature", negate = TRUE))
  
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
    
    # Calculate mean, standard deviation, and coeefficient of variation for groups of sequins
    # Create a logical vector determining if the sequin is within the threshold
    mutate(
      grouped_seq_cov = map(cov_tab, ~inner_join(., sequins, by = "Feature") %>%
                              select(Feature, Coverage, Concentration)),
      grouped_seq_cov = map2(grouped_seq_cov, Dilution,
                             ~mutate(.x, Concentration = Concentration/.y) %>%
                               group_by(Concentration) %>%
                               summarise(mean_cov = mean(Coverage),
                                         sd_cov = sd(Coverage)) %>%
                               na_if(0) %>% 
                               mutate(coe_var = sd_cov*100/mean_cov) %>%
                               mutate(threshold_detection = coe_var <= coe_of_variation))) %>%
    
    #Create a list of samples in which sequins were not detected
    mutate(under_detected = map(grouped_seq_cov, ~.x %>% filter(is.na(mean_cov)) %>% select(Concentration))
           ) %>%
  
    # perform linear regression on coverage vs conc., extract lm params, make plots
    mutate(
      seq_cov_filt = map2(seq_cov,grouped_seq_cov, ~ inner_join(.x, .y , by = "Concentration") %>%
                           filter(., Coverage > 0)),
      seq_cov_filt = map2(seq_cov_filt, lod, ~.x %>%
                            filter(Concentration >= .y) %>%
                            filter(., coe_var <= coe_of_variation) %>% #remove zero coverage values before lm
                            mutate(
                            lod = .y))) %>% 
    
    mutate(
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
                           geom_point(aes(shape = threshold_detection)) + 
                           geom_smooth(method = "lm") + 
                           stat_regline_equation(label.x= -0.1, label.y = 3) + 
                           stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = -0.1, label.y = 3.5) + 
                           xlab("Coverage (log[read depth])") + 
                           ylab("DNA Concentration (log[attamoles/uL])") + 
                           scale_shape(name = "Coefficient of variation", labels = c("below the threshold", "above the threshold")) +
                           theme_bw() 
                     ),
                     map(seq_cov_filt, # non-scaled plot if true
                         ~ ggplot(data=. , aes(x=Coverage, y= Concentration)) + 
                           geom_point(aes(color = threshold_detection)) + 
                           geom_smooth(method = "lm") + 
                           stat_regline_equation(label.x= 0, label.y = 1000) + 
                           stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 0, label.y = 100000) + 
                           xlab("Coverage (read depth)") + 
                           ylab("DNA Concentration (attamoles/uL)") +
                           scale_color_discrete(name = "Threshold detection") +
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
      mag_det = mag_ab,
      mag_ab = map(mag_ab, ~ select(., Feature, Concentration)), # drop Coverage column
      mag_ab = map2(mag_ab, Sample, ~ setNames(.x, c("Feature", .y))), #put sample name in MAG table
      
      # Remove MAGs below LOD
      #mag_det = mag_ab,
      #mag_det = map2(mag_det, Sample, ~ setNames(.x, .y, 'Concentration')), #get sample name out of header for filter
      mag_det = map2(mag_det, lod, ~ filter(.x, Concentration > .y)),
      mag_det = map(mag_det, ~select(.x, Feature, Concentration)),
      mag_det = map2(mag_det, Sample, ~ setNames(.x, c("Feature", .y))) #change header back to sample
    )

  # compile feature abundance across samples
  mag_tab <- scale_fac$mag_ab %>% 
    reduce(left_join, by="Feature") %>%
    column_to_rownames(var = "Feature")
  
  # extract plots
  plots <- scale_fac %>% select(Sample, plots)
  #Save scaling plots in .pdf format in the working directory
  plots <- plots %>%
    mutate(save_plots = map2(plots, Sample,  ~ggsave(filename = paste("",.y, ".pdf", sep=""), plot = .x, path = "../mock_output_data/")))
    
  
  # extract feature detection 
  mag_det <- scale_fac$mag_det %>% 
    reduce(left_join, by="Feature") %>%
    column_to_rownames(var = "Feature")
  
  
  # make list of results
  results <- list("mag_tab" = mag_tab, 
                  "mag_det" = mag_det, 
                  "plots" = plots, 
                  "scale_fac" = scale_fac)
  
  # return results
  attach(results)
  return(results)
  
}

