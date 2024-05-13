# Load required libraries
library(data.table)


# Command-line argument parsing
args <- commandArgs(trailingOnly = TRUE)
input_ncbiTax_path <- args[1]
input_ncbiComp_path <- args[2]
pubchemID <- args[3]


#Column 1: bsid of organism specific biosystem  
#Column 2: Taxonomy ID (TaxID) of source organism   
#Column 3: score 
dt_t = fread(input_ncbiTax_path)
colnames(dt_t) = c('bsid', 'TaxID', 'score')
dt_t[, bsid := as.integer(bsid)]

print(dt_t)

#Column 1: bsid of biosystem  
#Column 2: CID (compound ID) of the small molecule  
#Column 3: score 
dt_c = fread(input_ncbiComp_path)
colnames(dt_c) = c('bsid', 'cid', 'score')
dt_c[, bsid := as.integer(bsid)]
dt_c[, cid := as.integer(cid)]

print(dt_c)

print(pubchemID)

dt_ncbi_records = data.table(ID = unique(dt_t[bsid %in% unique(dt_c[cid == as.integer(pubchemID)]$bsid)]$TaxID))

print(dt_ncbi_records)

# Write results
fwrite(dt_ncbi_records, 'ncbi_records.csv')
