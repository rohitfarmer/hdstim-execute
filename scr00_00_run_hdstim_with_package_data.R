#!/usr/bin/R

# Purpose: To run HDStIM on CyTOF data included in the package.

library(HDStIM)

# NOTE: Type chi11 to see the built in package data

# Run the main HDStIM function using the built in package data
mapped_data <- HDStIM(chi11$expr_data, chi11$state_markers,
                       chi11$cluster_col, chi11$stim_label,
                       chi11$unstim_label, seed_val = 123, umap = TRUE,
                       umap_cells = 50, verbose = TRUE)


# Generate diagnostic plots
## These functions return ggplot objects. However, they can also generate figures
## if the path is given (not NULL).

cytof_k_plots <- plot_K_Fisher(mapped_data, path = NULL, verbose = TRUE)

cytof_u_plots <- plot_umap(mapped_data, path = NULL, verbose = TRUE)

cytof_e_plots <- plot_exprs(mapped_data, path = NULL, verbose = TRUE)

# Run marker ranking
## Marker ranking function returns a list with attribute statistics along with ggplot
## objects for figures. Similar to plotting functions above, if a path is given
## then the function can generate the plots as well.

marker_ranking <- marker_ranking_boruta(mapped_data,
  path = NULL,
  n_cells = NULL,
  max_runs = 100,
  seed_val = 123,
  verbose = 1
)

# Generate marker ranking heatmap
## This function will return a ComplexHeatmap object which can be printed
## on a chosen graphics device. 
pht <- plot_marker_ranking_heatmap(marker_ranking)


