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

# Read and process data
dt_redu <- fread(input_redu_path)
dt_redu <- dt_redu[grepl('|', NCBITaxonomy, fixed = TRUE)]
dt_redu[, NCBI := strsplit(unique(NCBITaxonomy), '|', fixed = TRUE)[[1]][1], by = .(NCBITaxonomy)]
dt_redu[, tax_name := strsplit(unique(NCBITaxonomy), '|', fixed = TRUE)[[1]][2], by = .(NCBITaxonomy)]
dt_redu[, NCBI := as.character(NCBI)]
dt_redu <- dt_redu[!is.na(NCBI) & NCBI != "" & !is.na(tax_name) & tax_name != ""]
dt_redu[, match_col := USI2MASST_matchCol(unique(USI)), by = .(USI)]
dt_redu <- dt_redu[, !"USI"]


fwrite(dt_redu, 'redu.tsv')