#!/usr/bin/R

# Purpose: To run marker ranking by Boruta.

library(tidyverse)
library(ComplexHeatmap)

results_folder <- file.path("results", "cytof-pilot12")
figures_folder <- file.path("figures", "cytof-pilot12")

cytof_m_ranks <- readRDS(file.path(results_folder, "marker-ranking.rds"))

cytof_mat <- cytof_m_ranks$attribute_stats %>%
  dplyr::select(stim_type, cell_population, state_marker, meanImp) %>%
  dplyr::group_by(stim_type, cell_population) %>%
  mutate(min_max = (meanImp - min(meanImp)) /(max(meanImp)-min(meanImp))) %>%
  dplyr::select(stim_type, cell_population, state_marker, min_max) %>%
  dplyr::mutate("stim_pop" = paste0(stim_type, "_", cell_population)) %>%
  dplyr::ungroup()

mat <- matrix(nrow = length(unique(cytof_mat$stim_pop)), ncol = length(unique(cytof_mat$state_marker)))
colnames(mat) <- unique(cytof_mat$state_marker)
rownames(mat) <- unique(cytof_mat$stim_pop)
for(i in 1:nrow(cytof_mat)){
        rowi <- cytof_mat$stim_pop[i]
        coli <- as.character(cytof_mat$state_marker[i])
        val <- cytof_mat$min_max[i]
        mat[rowi, coli] <- val 
}

library(circlize)
col_fun = colorRamp2(c(0, 1), c("black", "yellow"))
col_fun(seq(0, 1))

r_breaks <- c("LPS", "LPS", rep("PMA",9), rep("IFA", 9), rep("TCR" ,3), rep("IFG", 4))

hmap <- ComplexHeatmap::Heatmap(mat, cluster_rows = FALSE, cluster_columns = TRUE, 
                                col = col_fun, row_split = r_breaks,
                                row_title = NULL,
                                heatmap_legend_param = list(title = "Mean\nImportance"),
                                column_names_gp = gpar(fontsize = 10),
                                row_names_gp = gpar(fontsize = 10))
png(filename = file.path(figures_folder, "marker_ranking_heatmap.png") , width = 6, 
    height = 6, units = "in", pointsize = 12, bg = "white", res = 600)
draw(hmap)
dev.off()

# Z-score heatmaps.

cytof_mat_z <- cytof_m_ranks$attribute_stats %>%
  dplyr::select(stim_type, cell_population, state_marker, meanImp) %>%
  dplyr::group_by(stim_type, cell_population) %>%
  mutate(z_score = (meanImp - mean(meanImp)) / sd(meanImp)) %>%
  dplyr::select(stim_type, cell_population, state_marker, z_score) %>%
  dplyr::mutate("stim_pop" = paste0(stim_type, "_", cell_population)) %>%
  dplyr::ungroup()

mat_z <- matrix(nrow = length(unique(cytof_mat_z$stim_pop)), ncol = length(unique(cytof_mat_z$state_marker)))
colnames(mat_z) <- unique(cytof_mat_z$state_marker)
rownames(mat_z) <- unique(cytof_mat_z$stim_pop)
for(i in 1:nrow(cytof_mat_z)){
  rowi <- cytof_mat_z$stim_pop[i]
  coli <- as.character(cytof_mat_z$state_marker[i])
  val <- cytof_mat_z$z_score[i]
  mat_z[rowi, coli] <- val 
}


library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_fun(seq(-3, 0, 3))


hmap_z <- ComplexHeatmap::Heatmap(mat_z, cluster_rows = FALSE, cluster_columns = TRUE, 
                                col = col_fun, row_split = r_breaks,
                                row_title = NULL,
                                heatmap_legend_param = list(title = "Mean\nImportance\nZ-Score"),
                                column_names_gp = gpar(fontsize = 10),
                                row_names_gp = gpar(fontsize = 10))
png(filename = file.path(figures_folder, "marker_ranking_heatmap_z.png") , width = 6, 
    height = 6, units = "in", pointsize = 12, bg = "white", res = 600)
draw(hmap_z)
dev.off()
