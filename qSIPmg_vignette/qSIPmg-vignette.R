## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Load libraries-----------------------------------------------------------
#Load required libraries
library(tidyverse)
library(phyloseq)
library(HTSSIP)
library(ggpubr)
library(SIPmg)
set.seed(seed = 1000)

## ----Load data----------------------------------------------------------------
## Load data
#Coverage metadata
#Uncomment if your coverage data is in the format mentioned above for this file. Remains commented if you are using the output from `checkm coverage`
f_tibble <- readr::read_csv("mock_input_data/coverage_metadata.csv")

#Sequins metadata
sequins <- readr::read_csv(file="mock_input_data/sequins_metadata.csv")

#Dilutions data
seq_dil = readr::read_csv(file = "mock_input_data/dilutions_data.csv")

#Log scale BOOLEAN. True or False depending on how you would want the MAG coverages to be scaled. Select TRUE if you need MAG concentrations scaled on the log scale
log_scale = TRUE

#coe_of_variation. Acceptable coefficient of variation for coverage and detection (eg. 20 - for 20 % threshold of coefficient of variation) (Coverages above the threshold value will be flagged in the plots)
coe_of_variation = 20

#Taxonomy
gtdbtk_bac_summary = readr::read_delim("mock_input_data/gtdbtk.bac120.summary.tsv",
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
gtdbtk_archaea = readr::read_delim("mock_input_data/gtdbtk.ar122.summary.tsv",
                             "\t", escape_double = FALSE, trim_ws = TRUE)
#GC content
GC_content <- readr::read_csv(file = "mock_input_data/GC_content.csv")

#Fractions
fractions = readr::read_csv("mock_input_data/fractions.csv")


## ----estimate absolute concentrations-----------------------------------------
taxonomy_tibble = dplyr::bind_rows(gtdbtk_bac_summary, gtdbtk_archaea) #Combine bacteria and archaea taxonomy files if it has not been done yet
#mag_tab is a tibble with absolute concentrations of MAGs obtained by scaling MAG coverages using linear regression models on sequin coverages and concentration

##Scale MAG coverages to obtain MAG absolute concentrations and save scaling plots in the working directory
#For rlm scaling using scale_features_rlm
#For rlm scaling using scale_features_lm
mag_tab_scaled <- SIPmg::scale_features_rlm(f_tibble, sequins, seq_dil, log_scale, coe_of_variation = coe_of_variation)
mag_tab = as.matrix(mag_tab_scaled$mag_tab) #Extract absolute abundances as a matrix

## ----example scaling plot, echo = FALSE---------------------------------------
mag_tab_scaled$plots$plots[[1]]

## ----Make phyloseq objects----------------------------------------------------

mag.table = phyloseq::otu_table(mag_tab, taxa_are_rows = TRUE) #Phyloseq OTU table

taxonomy.object = SIPmg::tax.table(taxonomy_tibble) # Create a taxonomy phyloseq object
samples.object = SIPmg::sample.table(fractions) # Create a samples phyloseq object
phylo.qSIP = SIPmg::phylo.table(mag.table, taxonomy.object, samples.object) # Make a phyloseq table for downstream qSIP analysis

## ----Calculate atom fraction excess-------------------------------------------
atomX = SIPmg::qSIP_atom_excess_MAGs(phylo.qSIP,
                               control_expr='Isotope=="12C"',
                               treatment_rep='Replicate',
                               Gi = GC_content)
#Bootstrap confidence intervals
df_atomX_boot = SIPmg::qSIP_bootstrap_fcr(atomX, n_boot=10) #Change "parallel = FALSE" to compute using a single-core


## ----Plot atom fraction excess------------------------------------------------
CI_threshold = 0
df_atomX_boot = df_atomX_boot %>%
  dplyr::mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))

(atom_f_excess_plot = ggplot2::ggplot(df_atomX_boot, aes(OTU, A, ymin=A_CI_low, ymax=A_CI_high, color=Incorporator)) +
  geom_pointrange(size=0.25) +
  geom_linerange() +
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5) +
  labs(x='MAGs', y='Atom fraction excess') +
  theme_bw() +
  coord_flip() +
  ggtitle("Isotope incorporating MAGs"))
ggplot2::ggsave(filename = "atom_fration_excess.pdf", plot = atom_f_excess_plot)

## ----incorporators------------------------------------------------------------

#Get incorporator info
n_incorp = df_atomX_boot %>%
  dplyr::filter(Incorporator == TRUE) %>%
  nrow
#Get incorporator list
incorporator_list = SIPmg::incorporators_taxonomy(taxonomy = taxonomy_tibble, bootstrapped_AFE_table = df_atomX_boot)
#Print incorporator information
cat('Number of incorporators:', n_incorp, '\n')
print(incorporator_list, n = nrow(incorporator_list))

## ----SIP plots,comment=FALSE, message=FALSE, warning=FALSE--------------------
#Load function for abundance plots
data_vis_table = SIPmg::sip_vis(fractions = fractions, mag_tab = mag_tab, afe_tab = df_atomX_boot, taxonomy = taxonomy_tibble)

(plot_bd = ggplot2::ggplot(data_vis_table, aes(x = buoyant_density, y = normalized_abundance,
                        shape = replicate, color = isotope)) +
  geom_point() +
  #scale_alpha_discrete(range = c(1,0.4)) +
  scale_color_manual(values=c("gold", "deepskyblue2")) +
  geom_line() +
  theme_bw() +
  xlab("Buoyant density (g/mL)") +
  ylab("Normalized abundance (attaM)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #scale_x_continuous(breaks = x_ticks, limits = c(1.69, 1.73)) +
  facet_wrap(~taxonomy, ncol = 4))

# Comparison of centroid buoyant densities of labeled and unlabeled MAGs

(ggplot2::ggplot(data_vis_table %>%
  dplyr::group_by(feature, isotope) %>%
  dplyr::mutate(centroid_BD = stats::weighted.mean(x = buoyant_density,
                                     w = absolute_abundance)) , aes(isotope, centroid_BD, fill = isotope)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~taxonomy, ncol = 4))


## ----Load data 2, echo = FALSE------------------------------------------------
##Load data
#Coverage data
f_tibble <- readr::read_csv("mock_input_data/coverages_outliers.csv")

#Sequins metadata
sequins <- readr::read_csv(file="mock_input_data/sequin_meta_outliers.csv")

#Dilutions data
seq_dil = readr::read_csv(file = "mock_input_data/seq_dilution_outliers.csv")

#Log scale BOOLEAN. True or False depending on how you would want the MAG coverages to be scaled. Select TRUE if you need MAG concentrations scaled on the log scale
log_scale = TRUE

#coe_of_variation. Acceptable coefficient of variation for coverage and detection (eg. 20 - for 20 % threshold of coefficient of variation) (Coverages above the threshold value will be flagged in the plots)
coe_of_variation = 50


## ----robust linear regression-------------------------------------------------
mag_tab_scaled_rlm <- SIPmg::scale_features_rlm(f_tibble, sequins, seq_dil, log_scale, coe_of_variation = coe_of_variation)

## ----linear regression--------------------------------------------------------
mag_tab_scaled_lm <- SIPmg::scale_features_lm(f_tibble, sequins, seq_dil, log_scale, coe_of_variation = coe_of_variation, cook_filtering = TRUE)

## ----load images, echo = F----------------------------------------------------
rlm_example = EBImage::readImage("rlm-example.png")
rlm_example = EBImage::resize(rlm_example,dim(rlm_example)[1]/2)

lm_example = EBImage::readImage("lm-example.png")
lm_example = EBImage::resize(lm_example,dim(lm_example)[1]/2)

filtered_lm_example = EBImage::readImage("filtered_lm-example.png")
fitlered_lm_example = EBImage::resize(filtered_lm_example,dim(filtered_lm_example)[1]/2)

cooksd_example = EBImage::readImage("cooksd-example.png")
cooksd_example = EBImage::resize(cooksd_example,dim(cooksd_example)[1]/2)


## ----show images--------------------------------------------------------------
# Robust linear regression plot
EBImage::display(rlm_example)

#Linear regression plot without filtering sequin data
EBImage::display(lm_example)

## ---- outlier sequin coverages------------------------------------------------

#Cook's distance threshold of the data set
4/(length(mag_tab_scaled_lm$scale_fac$cooksd[[3]]))

#Cook's distance of the outlier before filtering the data
max(mag_tab_scaled_lm$scale_fac$cooksd[[3]])

#Before filtration
EBImage::display(cooksd_example)

## ----better fit---------------------------------------------------------------
#Linear regression plot with filtered outliers in sequin data
EBImage::display(filtered_lm_example)

## ----Load data 3, echo = F----------------------------------------------------
## Load data
#Coverage metadata
#Uncomment if your coverage data is in the format mentioned above for this file. Remains commented if you are using the output from `checkm coverage`
f_tibble <- readr::read_csv("mock_input_data/coverage_metadata.csv")

#Sequins metadata
sequins <- readr::read_csv(file="mock_input_data/sequins_metadata.csv")

#Dilutions data
seq_dil = readr::read_csv(file = "mock_input_data/dilutions_data.csv")

#Log scale BOOLEAN. True or False depending on how you would want the MAG coverages to be scaled. Select TRUE if you need MAG concentrations scaled on the log scale
log_scale = TRUE

#coe_of_variation. Acceptable coefficient of variation for coverage and detection (eg. 20 - for 20 % threshold of coefficient of variation) (Coverages above the threshold value will be flagged in the plots)
coe_of_variation = 20

#Taxonomy
gtdbtk_bac_summary = readr::read_delim("mock_input_data/gtdbtk.bac120.summary.tsv",
                                 "\t", escape_double = FALSE, trim_ws = TRUE)
gtdbtk_archaea = readr::read_delim("mock_input_data/gtdbtk.ar122.summary.tsv",
                             "\t", escape_double = FALSE, trim_ws = TRUE)
#GC content
GC_content <- readr::read_csv(file = "mock_input_data/GC_content.csv")

#Fractions
fractions = readr::read_csv("mock_input_data/fractions.csv")


## ----estimate absolute concentrations 2, echo = F-----------------------------
taxonomy_tibble = dplyr::bind_rows(gtdbtk_bac_summary, gtdbtk_archaea) #Combine bacteria and archaea taxonomy files if it has not been done yet
#mag_tab is a tibble with absolute concentrations of MAGs obtained by scaling MAG coverages using linear regression models on sequin coverages and concentration

##Scale MAG coverages to obtain MAG absolute concentrations and save scaling plots in the working directory
#For rlm scaling using scale_features_rlm
#For rlm scaling using scale_features_lm
mag_tab_scaled <- SIPmg::scale_features_rlm(f_tibble, sequins, seq_dil, log_scale, coe_of_variation = coe_of_variation)
mag_tab = as.matrix(mag_tab_scaled$mag_tab) #Extract absolute abundances as a matrix

## ----Make phyloseq objects 2, echo = F----------------------------------------

mag.table = phyloseq::otu_table(mag_tab, taxa_are_rows = TRUE) #Phyloseq OTU table

taxonomy.object = SIPmg::tax.table(taxonomy_tibble) # Create a taxonomy phyloseq object
samples.object = SIPmg::sample.table(fractions) # Create a samples phyloseq object
phylo.qSIP = SIPmg::phylo.table(mag.table, taxonomy.object, samples.object) # Make a phyloseq table for downstream qSIP analysis

## ----AFE methodGet bootstrapped AFE table 2-----------------------------------
atomX = SIPmg::qSIP_atom_excess_MAGs(phylo.qSIP,
                               control_expr='Isotope=="12C"',
                               treatment_rep='Replicate',
                               Gi = GC_content)
#Bootstrap confidence intervals
df_atomX_boot = SIPmg::qSIP_bootstrap_fcr(atomX, n_boot=10) #Change "parallel = FALSE" to compute using a single-core
CI_threshold = 0
df_atomX_boot = df_atomX_boot %>%
  dplyr::mutate(Incorporator = A_CI_low > CI_threshold,
         OTU = reorder(OTU, -A))

## ----list incorporators 2-----------------------------------------------------
#Get incorporator info
n_incorp_AFE = df_atomX_boot %>%
  dplyr::filter(Incorporator == TRUE) %>%
  nrow
#Get incorporator list
incorporator_list_AFE = SIPmg::incorporators_taxonomy(taxonomy = taxonomy_tibble, bootstrapped_AFE_table = df_atomX_boot)
#Print incorporator information
cat('Number of incorporators:', n_incorp_AFE, '\n')
rmarkdown::paged_table(incorporator_list_AFE)

## -----------------------------------------------------------------------------
sessionInfo()

