#!/usr/bin/R

# Purpose: To run HDStIM on CyTOF data.

library(tidyverse)
library(HDStIM)
library(arrow)

results_folder <- file.path("results", "cytof-pilot12")
figures_folder <- file.path("figures", "cytof-pilot12")
dir.create(figures_folder, recursive = TRUE)


# Load data.
cytof_dat <- arrow::read_feather(file.path(results_folder,"cytof-dat-asinh-transform.feather"))

# Read pannel info
panel <- read_tsv(file.path("meta", "panel", "cytof_phospho_panel_v4.txt"), show_col_types = FALSE)

stims <- setdiff(unique(cytof_dat$stim_type), "UNS")
unstim <- "UNS"
state_marker <- panel %>% dplyr::filter(marker_class == "state") %>%
        dplyr::pull(antigen)
state_marker <- str_replace_all(state_marker, "-", "_")

cytof_mapped_data <-  HDStIM(cytof_dat, state_marker,
                             "cell_population", stims, unstim, seed_val = 123, 
                             umap = TRUE, umap_cells = 5000, 
                             verbose = TRUE)

saveRDS(cytof_mapped_data, file.path(results_folder, "cytof-mapped-data.rds"))

# Generate diagnostic plots
cytof_k_plots <- plot_K_Fisher(cytof_mapped_data, path = file.path(figures_folder, "cytof-k-Fisher-plots") , verbose = TRUE)

cytof_u_plots <- plot_umap(cytof_mapped_data, path = file.path(figures_folder, "cytof-u-plots"), verbose = TRUE)

cytof_e_plots <- plot_exprs(cytof_mapped_data, path = file.path(figures_folder, "cytof-exprs-plots"), verbose = TRUE)

# Run marker ranking
marker_ranking <- marker_ranking_boruta(cytof_mapped_data,
  path = figures_folder,
  n_cells = NULL,
  max_runs = 100,
  seed_val = 123,
  verbose = 1
)
saveRDS(ranking, file.path(results_folder, "marker-ranking.rds"))

# Generate marker ranking heatmap
pht <- plot_marker_ranking_heatmap(marker_ranking)

