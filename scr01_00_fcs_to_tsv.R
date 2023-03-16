#!/usr/bin/R

# Purpose: To read FCS files and export expression data with additional columns 
# to tsv format to be used with HDStIM

library(flowCore)
library(tidyverse)
library(arrow)

data_folder <- file.path("data", "cytof-pilot12")
results_folder <- file.path("results", "cytof-pilot12")
dir.create(results_folder,recursive = TRUE )

# CYTOF
cytof_files_fcs <- read_tsv(file.path("data", "cytof-pilot12","cytof-pilot12.tsv"), 
                            show_col_types = FALSE)
panel <- read_tsv(file.path("meta", "panel", "cytof_phospho_panel_v4.txt"), show_col_types = FALSE)
cytof_dat_out <- tibble()
for(j in 1:nrow(cytof_files_fcs)){
        cytof_file_path <- file.path(data_folder, cytof_files_fcs$fcs_file[j])
        cytof_fcs_dat <- read.FCS(cytof_file_path, transformation = FALSE, truncate_max_range = FALSE, alter.names=TRUE)
        cytof_expr_dat <- exprs(cytof_fcs_dat) %>%
                as_tibble()
        cytof_col_names <- as.character(cytof_fcs_dat@parameters@data$name) # Channel names in the FCS columns.
        cytof_col_desc <- as.character(cytof_fcs_dat@parameters@data$desc) # Name of the markers. 
        c_names <-  coalesce(cytof_col_desc, cytof_col_names) %>% # User marker names where available else use channel names for example,
        str_split( "_", simplify = TRUE)                          # for time and forward and side scattering.
        c_names <- gsub("-", "_", coalesce(na_if(c_names[,2], ""), c_names[,1]))
        features <- c(gsub("-", "_", panel$antigen), "Time")
        colnames(cytof_expr_dat) <- c_names
        cytof_expr_dat <- dplyr::select(cytof_expr_dat, all_of(features))

        cytof_dat_out <- bind_rows(cytof_dat_out, tibble("cell_population" = cytof_files_fcs$cell_population[j], 
                                                         "stim_type" = cytof_files_fcs$stim_type[j], 
                                                         "condition" = cytof_files_fcs$condition[j],
                                                         "pilot" = cytof_files_fcs$pilot[j],
                                                         cytof_expr_dat))
}
#arrow::write_feather(cytof_dat_out, file.path(results_folder, "cytof-dat-no-transform.feather"))

# Arch sinh transform
markers <- str_replace(panel$antigen, "-", "_")

asi <- function(x) asinh(x/5)

cytof_dat_out[markers] <- apply(cytof_dat_out[markers],2, asi)


arrow::write_feather(cytof_dat_out, file.path(results_folder, "cytof-dat-asinh-transform.feather"))
