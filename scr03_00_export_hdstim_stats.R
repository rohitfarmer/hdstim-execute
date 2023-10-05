#!/usr/bin/R

# Purpose: Export stats from the HDStIM output.

library(tidyverse)

results_folder <- file.path("results", "all-pilots-subsets")
figures_folder <- file.path("figures", "all-pilots-subsets")

cytof_mapped_data <- readRDS(file.path(results_folder, "cytof-mapped-data.rds"))

# Fisher's exact test results.
write_tsv(cytof_mapped_data$all_fisher_p_val, file.path(results_folder, "cytof-fishers-data.tsv"))



