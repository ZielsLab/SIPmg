library(tidyverse)
#-----------------------------------------------------------------------
#' This function provides the plots for abundance vs Buoyant density of incorporators.
#' The plots are saved in the directory "abundance_plots" in the current working directory 
#' of the "qSIP_analysis.rmd" pipeline
#' 
#' Input parameters:
#' - incorporator_list - A tibble consisting of two columns. This tibble is 
#' obtained from the `qSIP_analysis.rmd` R markdown
#' First column consists of incorporator feature ID.
#' Second column consists of the incorporator taxonomy
#' 
#' - fractions: Fractions metadata file loaded in as part of the `qSIP_analysis.rmd` R markdown
#' 
#' - mag_tab: A tibble of feature abundances obtained from the `qSIP_analysis.rmd` R markdown
#' Each column consists of the scaled abundances of the features. Row names idenitfy the Feature
#' 
#' Output:
#'  - "mean_abundance_" prefix is for the plots visualizing the abundance estimates of features
#'  over the buoyant density gradient
#'  - "Rep_1_abundance" prefix is for the plots visualizing the abundance estimates of features
#'  from replicate 1 over the buoyant density gradient
#'  - "Rep_2_abundance" prefix is for the plots visualizing the abundance estimates of features
#'  from replicate 2 over the buoyant density gradient
#'  - "Rep_3_abundance" prefix is for the plots visualizing the abundance estimates of features
#'  from replicate 3 over the buoyant density gradient
#'  
#'  The goal of plotting individual replicates separately is to identify any anomalies,
#'  if any, regarding the changes in the scaled abundances across the buoyant density gradient
#'  which may not be obvious in the plots with the mean and standard deviation of abundances. 
#'  These plots could help in debugging
#------------------------------------------------------------------------

plot_abundance = function(incorporator_list, fractions, mag_tab) {
dir.create("abundance_plots")
#Rename first column of `incorporator_list_rename` tibble to Feature from OTU
incorporator_list_renamed = incorporator_list %>%
  rename("Feature" = "OTU")
#Reformat mag_tab to make rownames in this tibble as column names with the column name Feature
mag_tab_reformatted = as.data.frame(mag_tab) %>%
  rownames_to_column(var = "Feature")
#Make a tibble by joining abundance and MAG names from incorporator list
incorporator_abs_abundance= incorporator_list_renamed %>%
  inner_join(., mag_tab_reformatted, by = "Feature")
#Pivot the tibble to later access abundance values of every Bin in every sample
#Store this as a temporary tibble
long_cov = incorporator_abs_abundance %>%
  pivot_longer(fractions$Sample, names_to = "Sample", values_to = "Abundance")
#Keep Feature and taxonomy columns only
incorporator_abs_abundance = incorporator_abs_abundance %>%
  select(Feature, Taxonomy)

#Create a tibble within the previous tibble to access abundances of each MAG in every sample
incorporator_abs_abundance = incorporator_abs_abundance %>%
  mutate(
    abs_abundance = map(Feature, ~ filter(long_cov, Feature == .) %>% select(Sample, Abundance))
  )
#Populate the corresponding fraction, BD, Replicate, and Isotope information
#for all samples. This gives coverage values for a certain replicate, in a certain BD fraction, for a particular isotope treatment
#for all MAGs
incorporator_abs_abundance = incorporator_abs_abundance %>%
  mutate(
    abs_abundance = map(abs_abundance, ~ mutate(., Isotope = fractions$Isotope,
                                              Fraction = fractions$Fraction,
                                              Buoyant_density = fractions$Buoyant_density,
                                              Rep = str_split_fixed(fractions$Sample, pattern = "_", n= Inf)[,3]))
  )
#Summarise the mean abundance, standard deviation of abundance, and mean BD for each fraction and isotope treatment
incorporator_abs_abundance = incorporator_abs_abundance %>%
  mutate(
    summary_coverage = map(abs_abundance, ~ group_by(., Fraction, Isotope ) %>%
                           summarise(mean_abs_abundance = mean(Abundance),
                                     mean_BD = mean(Buoyant_density),
                                     sd_abs_abundance = sd(Abundance))),
    summary_coverage = map(summary_coverage, ~ arrange(., Isotope))
                             
  )
#Plot abundance vs BD and save the plots in the current path
incorporator_abs_abundance = incorporator_abs_abundance %>%
  mutate(
    plots = map(summary_coverage, ~ ggplot(data = ., aes(x = mean_BD, y = mean_abs_abundance)) +
                  geom_point(aes(color = Isotope)) +
                  geom_line(aes(color = Isotope)) +
                  geom_errorbar(aes(ymin = mean_abs_abundance - sd_abs_abundance,
                                    ymax = mean_abs_abundance + sd_abs_abundance)) +
                  ylab("Mean absolute \n abundance (attamole/uL)") +
                  xlab("Mean buoyant \n density (g/mL)") +
                  theme_bw()),
    save_plots = map2(plots, Feature,  ~ggsave(filename = paste("mean_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/"))
  )
#' View plots of abundance vs buoyant density for each replicate of the incorporators.
#' This could give a better visualization of anomalies, if any, present in the
#' abundance vs BD plots which may not be obvious in the mean and standard deviation estimates

incorporator_abs_abundance = incorporator_abs_abundance %>%
  mutate(
    plots_Rep1 = map(abs_abundance, ~ filter(., Rep == 1) %>%
                  ggplot(data = ., aes(x = Buoyant_density, y = Abundance)) +
                  geom_point(aes(color = Rep, shape = Isotope)) +
                  geom_line(aes(color = Isotope)) +
                  ylab("Absolute \n abundance (attamole/uL)") +
                  xlab("Buoyant \n density (g/mL)") +
                  theme_bw()),
    plots_Rep2 = map(abs_abundance, ~ filter(., Rep == 2) %>%
                       ggplot(data = ., aes(x = Buoyant_density, y = Abundance)) +
                       geom_point(aes(color = Rep, shape = Isotope)) +
                       geom_line(aes(color = Isotope)) +
                       ylab("Absolute \n abundance (attamole/uL)") +
                       xlab("Buoyant \n density (g/mL)") +
                       theme_bw()),
    plots_Rep3 = map(abs_abundance, ~ filter(., Rep == 3) %>%
                       ggplot(data = ., aes(x = Buoyant_density, y = Abundance)) +
                       geom_point(aes(color = Rep, shape = Isotope)) +
                       geom_line(aes(color = Isotope)) +
                       ylab("Absolute \n abundance (attamole/uL)") +
                       xlab("Buoyant \n density (g/mL)") +
                       theme_bw()),
    save_plots = map2(plots_Rep1, Feature,  ~ggsave(filename = paste("Rep_1_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/")),
    save_plots = map2(plots_Rep2, Feature,  ~ggsave(filename = paste("Rep_2_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/")),
    save_plots = map2(plots_Rep3, Feature,  ~ggsave(filename = paste("Rep_3_abundance_",.y, ".pdf", sep=""), plot = .x, path = "abundance_plots/"))
  )
}
