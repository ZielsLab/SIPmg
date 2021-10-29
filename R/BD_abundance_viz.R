#' Create a table for easy data visualization of plots with
#' buoyant densities on X-axis and absolute abundance on Y-axis
#'
#' @param fractions_df
#' A fractions file with the following columns
#' - Replicate: Depends on how many replicates the study has
#' - Fractions: Typically in the range of 2-24
#' - Buoyant_density: As calculated from the refractometer for each fraction and replicate
#' - Isotope: "12C", "13C", "14N", "15N" etc.
#' - DNA_concentration
#' - Sample: In the format "'isotope'_rep_#_fraction_#".
#'   For instance, "12C_rep_1_fraction_1"
#' @param mag_tab
#' mag_tab: a tibble with first column "Feature" that contains the bin ID (or contig IDs),
#' and the rest of the columns represent samples with features' scaled abundances (attamoles/uL)
#' @import magrittr
#' @importFrom rlang .data
#' @return  data  frame: formatted table for BD vs absolute abundance visualization
#' @export

sip_vis = function(fractions, mag_tab, taxonomy, afe_tab)
{

  #Retain incorporator IDs in a separate dataframe
  incorporator_list = incorporators_taxonomy(taxonomy, df_atomX_boot) %>%
    dplyr::rename("feature" = "OTU")

  #Make a data visualization table that has Substrate, replicate,
  #fraction,  buoyant density, absolute abundance,
  #and absolute abundance normalized to peak in the sample

  data_vis_table =  as.data.frame(mag_tab) %>%
    tibble::rownames_to_column(var = "feature") %>%
    dplyr::inner_join(incorporator_list) %>%
    tidyr::pivot_longer(-c(feature,taxonomy), names_to = "Sample", values_to = "absolute_abundance") %>%
    dplyr::inner_join(fractions) %>%
    tidyr::separate(Sample, into = c("isotope", "replicate", "rep_num", "fraction", "fraction_num"), sep = "_") %>%
    tidyr::unite("replicate", replicate:rep_num) %>%
    tidyr::unite("fraction", fraction:fraction_num) %>%
    dplyr::group_by(isotope, replicate, feature) %>%
    dplyr::rename("buoyant_density" = "Buoyant_density") %>%
    dplyr::mutate(normalized_abundance = absolute_abundance/max(absolute_abundance))

  return(data_vis_table)
}
