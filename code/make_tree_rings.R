# Load required libraries
library(data.table)
library(tools)

# Define function to generate match columns
USI2MASST_matchCol <- function(x){
    elements = strsplit(x, ':', fixed = TRUE)[[1]]
    mzspec = elements[1]
    msvid = elements[2]
    filename = tools::file_path_sans_ext(basename(elements[3]))
    paste0(c(mzspec, msvid, filename), collapse = ':')
}

# Command-line argument parsing
args <- commandArgs(trailingOnly = TRUE)
input_redu_path <- args[1]
input_masst_path <- args[2]
output_path <- args[3]

# Read and process data
dt_redu <- fread(input_redu_path)
dt_redu <- dt_redu[grepl('|', NCBITaxonomy, fixed = TRUE)]
dt_redu[, ncbiid := strsplit(unique(NCBITaxonomy), '|', fixed = TRUE)[[1]][1], by = .(NCBITaxonomy)]
dt_redu[, ID := as.character(ncbiid)]
dt_redu <- dt_redu[!is.na(ncbiid)]
dt_redu[, match_col := USI2MASST_matchCol(unique(USI)), by = .(USI)]
dt_redu <- dt_redu[, !"USI"]

dt_masst_results <- fread(input_masst_path)
dt_masst_results[, match_col := paste0(strsplit(unique(USI), ':', fixed = TRUE)[[1]][1:3], collapse = ':'), by = .(USI)]
dt_masst_results <- dt_masst_results[, .(Cosine = max(Cosine, na.rm = TRUE), `Matching Peaks` = `Matching Peaks`[which.max(Cosine)]), by = .(match_col)]

dt_masst_results <- dt_masst_results[dt_redu, on = .(match_col)]
dt_masst_results <- dt_masst_results[, .(n_detected = sum(!is.na(Cosine)), n_samples = .N, p_detected = sum(!is.na(Cosine))/.N*100), by = .(ncbiid)]
dt_masst_results[, ID := as.character(ncbiid)]

# Write results
fwrite(dt_masst_results, output_path)
