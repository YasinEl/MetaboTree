#!/usr/bin/env Rscript

# Load necessary libraries
library(tidyverse)
library(taxize)
library(data.table)



# Function to fetch classification with retries
get_classification_with_retries <- function(ncbi_id, max_attempts = 5) {
  for (attempt in 1:max_attempts) {
    classification_result <- tryCatch({
      #Sys.sleep(0.05)
      classification(ncbi_id, db = 'ncbi', batch_size = 4)
    }, error = function(e) {
      message(paste("Attempt", attempt, "failed for NCBI ID", ncbi_id, ": ", e$message))
      if (attempt < max_attempts) {
        Sys.sleep(2)  # Wait for 2 seconds before retrying
      }
      NULL
    })
    if (!is.null(classification_result)) {
      message(paste("Successfully fetched data on attempt", attempt, "for NCBI ID", ncbi_id))
      return(classification_result)
    }
  }
  message(paste("All attempts failed for NCBI ID", ncbi_id))
  return(NULL)
}



# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
redu_ncbi <- args[2]
output_file <- args[3]

# Set global options
Sys.setenv(ENTREZ_KEY = "db21f3d8fcc51ae12b1b0aef7710673c0508")
taxize_options(ncbi_sleep = 0.102, taxon_state_messages = TRUE, quiet = FALSE)


data_redu_linages <- fread(redu_ncbi)

# Read data from file
data <- fread(input_file)
data[, ncbiid := str_extract(NCBITaxonomy, "^[^|]*")]

# Filter and prepare NCBI IDs
ncbiid <- na.omit(unique(data[!(ncbiid %in% unique(data_redu_linages$ncbi))]$ncbiid))


if (length(ncbiid) > 0){


ncbiid <- as.numeric(ncbiid)
ncbiid <- ncbiid[!is.na(ncbiid)]


# Initialize variables for loop
collect_dt <- list()
requested_ids <- character()

# Main loop to fetch data
for (i in seq_along(ncbiid)) {
  current_id <- as.character(ncbiid[i])
  if (!current_id %in% requested_ids) {
    print(current_id)
    dt_li <- get_classification_with_retries(ncbiid[i])
    if (!is.null(dt_li)) {
      dt <- as.data.table(dt_li[[1]])
      dt[, ncbi := current_id]
      print(dt)
      collect_dt[[length(collect_dt) + 1]] <- dt
      requested_ids <- c(requested_ids, current_id)
    }
  } else {
    message(paste("Skipping already fetched entry at index", i))
  }
  print(paste(i, '/', length(ncbiid)))
}

# Combine results and write to file
dt_redu_linages <- rbindlist(collect_dt)

data_redu_linages = rbindlist(list(dt_redu_linages, data_redu_linages))

}


# Compute frequency of each name within each rank
data_redu_linages[, freq := .N, by = .(rank, name)]

# Sort by frequency and then by name to ensure selection of the less common name
setorder(data_redu_linages, rank, freq, name)

# Remove duplicate ranks within each ncbi group, keeping the less common name
data_redu_linages <- data_redu_linages[, .SD[1], by = .(ncbi, rank)]

# Transform the data.table using dcast
data_redu_linages <- dcast(data_redu_linages, ncbi ~ rank, value.var = "name")

setnames(data_redu_linages, 'ncbi', 'NCBI')

fwrite(data_redu_linages, output_file)