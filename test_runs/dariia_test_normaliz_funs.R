source("normaliz_funs.R")
# import mock data
cov_tibble <- read_csv("./mock_import/parsed_raw_cov.csv")
bin_tibble <- read_csv("./mock_import/fake_parsed_bin.csv")
opr <- "mean"

# test pooling function
pooled_bins <- pool_bin_stat(cov_tibble, bin_tibble, operator = opr)
write.csv(pooled_bins, "./mock_export/pool_bin_stat.csv", row.names = F)
