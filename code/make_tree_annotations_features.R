# Load required libraries
library(optparse)
library(ggtree)
library(data.table)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
library(ape)
library(ggnewscale)
library(grid)
library(magick)
library(stringr)
library(cowplot)
library(png)
library(grid)
library(readr)
library(viridis)

USI2MASST_matchCol <- function(x) {
  elements <- strsplit(x, ":")[[1]]
  mzspec <- elements[1]
  msvid <- elements[2]
  filename <- tools::file_path_sans_ext(basename(elements[3]))
  return(paste(mzspec, msvid, filename, sep = ":"))
}


# Setup command-line options
option_list <- list(
  make_option(c("-i", "--input_tree"), type="character", default=NULL, help="Input Newick file path", metavar="character"),
  make_option(c("-l", "--input_masst"), type="character", default=NULL, help="Input library file path", metavar="character"),
  make_option(c("-r", "--input_redu"), type="character", default=NULL, help="Input REDU file path", metavar="character"),
  make_option(c("-n", "--input_linage"), type="character", default=NULL, help="Input linage file path", metavar="character"),
  make_option(c("-u", "--usi"), type="character", default=NULL, help="USI", metavar="character"),
  make_option(c("-p", "--cid"), type="character", default=NULL, help="CID", metavar="character")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

lib_id = strsplit(args$usi, ':')[[1]][5]
cid = args$cid

# Load ReDU table
dt_redu <- fread(args$input_redu)

# get all columnnames sparql_* 
sparql_columns = colnames(dt_redu)[grep('^sparql_', colnames(dt_redu))]

dt_redu_match = dt_redu[!is.na(uid_leaf), c(.(NCBIDivision = NCBIDivision[1], UBERONBodyPartName = paste0(sort(unique(UBERONBodyPartName)), collapse= ', '), tax_name = tax_name[1]), 
  lapply(.SD, function(col) col[1])
), by =.(uid_leaf), .SDcols = sparql_columns]


# Load ReDU table
dt_linage <- fread(args$input_linage)
dt_linage = unique(dt_linage[, !c('uid')], by ='uid_leaf')
dt_redu_match = merge(dt_redu_match, dt_linage, by = 'uid_leaf', all.x = TRUE)

# Load tree
tree <- read.newick(args$input_tree)

# Get Tree ID base table
dt_tree_ids = data.table(uid_leaf = c(tree$tip.label, tree$node.label))

dt_tree_ids_redu <- merge(dt_tree_ids, dt_redu_match, by = "uid_leaf", all.x = TRUE)


# Load data tables
dt_masst <- fread(args$input_masst)


if(nrow(dt_masst) > 0){

  dt_masst = dt_masst[smiles_name != '']
  


  dt_masst[, match_col := USI2MASST_matchCol(USI[1]), by =.(USI)]


  max_indices <- dt_masst[, .I[which.max(Cosine)], by = .(match_col, smiles_name, scan_id)]$V1
  dt_masst <- dt_masst[max_indices]

  dt_masst <- dt_masst[, .(
    
    Cosine = max(Cosine), 
    `Matching Peaks` = max(`Matching Peaks`[which.max(Cosine)]), 
    feature_intensity_TICnorm_best = max(feature_intensity_TICnorm[which.max(Cosine)], na.rm = TRUE), 
    feature_intensity_CLR_best = max(feature_intensity_CLR[which.max(Cosine)], na.rm = TRUE), 
    feature_intensity_TICnorm_sum = sum(feature_intensity, na.rm= TRUE)/unique(feature_area_TIC), 
    feature_count = sum(feature_intensity_TICnorm[!is.na(feature_intensity_TICnorm) & feature_intensity_TICnorm > 0], na.rm = TRUE), 
    inchi_key_first_block_count_by_usi = paste0(unique(inchi_key_first_block), collapse=',')
    
    ), by = .(match_col, smiles_name)]

  dt_masst <- dcast(dt_masst, match_col ~ smiles_name, 
                   value.var = c(
                    "Cosine", "Matching Peaks", "inchi_key_first_block_count_by_usi",
                    "feature_intensity_TICnorm_best", 
                    "feature_intensity_CLR_best", 
                    "feature_intensity_TICnorm_sum", 
                    "feature_count"
                    ),)
  
  dt_redu_match = unique(dt_redu[, c('match_col', 'uid_leaf')])

  dt_redu_match <- merge(dt_masst, dt_redu_match, by = "match_col", all.x = TRUE, all.y = FALSE)

  #dt_redu_match = unique(dt_redu_match)

#print('unique2 done')

  dt_tree_ids_redu <- merge(dt_tree_ids_redu, dt_redu_match, by = "uid_leaf", all.x = TRUE, all.y = FALSE)

  dt_tree_ids_redu[, FeatureID := uid_leaf]

  #dt_tree_ids_redu = dt_tree_ids_redu[, .(Cosine = max(Cosine, na.rm = TRUE), `Matching Peaks` = `Matching Peaks`[which.max(`Matching Peaks`)], detected = any(detected), inchi_key_first_block_count = length(unique(unlist(unlist(strsplit(unique(inchi_key_first_block_count_by_usi), ',')))))), by =.(FeatureID, tax_name, UBERONBodyPartName, NCBIDivision, Database, biotype,clade,class,cohort,family,genus,infraclass,infraorder,isolate,kingdom,order,parvorder,phylum,section,species,`species group`,`species subgroup`,strain,subclass,subcohort,subfamily,subgenus,subkingdom,suborder,subphylum,subsection,subspecies,subtribe,superclass,superfamily,superkingdom,superorder,tribe,varietas)]

all_columns <- names(dt_tree_ids_redu)

print('all_columns_for_suffixes')
print(all_columns)

cosine_suffixes <- unique(sub("Cosine_", "", grep("^Cosine_", all_columns, value = TRUE)))
matching_peaks_suffixes <- unique(sub("Matching Peaks_", "", grep("^Matching Peaks_", all_columns, value = TRUE)))
inchi_key_suffixes <- unique(sub("inchi_key_first_block_count_by_usi_", "", grep("^inchi_key_first_block_count_by_usi_", all_columns, value = TRUE)))
feature_intensity_TICnorm_best_suffixes <- unique(sub("feature_intensity_TICnorm_best_", "", grep("^feature_intensity_TICnorm_best_", all_columns, value = TRUE)))
feature_intensity_CLR_best_suffixes <- unique(sub("feature_intensity_CLR_best_", "", grep("^feature_intensity_CLR_best_", all_columns, value = TRUE)))
feature_intensity_TICnorm_sum_suffixes <- unique(sub("feature_intensity_TICnorm_sum_", "", grep("^feature_intensity_TICnorm_sum_", all_columns, value = TRUE)))
feature_count_suffixes <- unique(sub("feature_count_", "", grep("^feature_count_", all_columns, value = TRUE)))



# Ensure we have a consistent list of suffixes across all patterns
suffixes <- intersect(intersect(cosine_suffixes, matching_peaks_suffixes), inchi_key_suffixes)
suffixes <- intersect(suffixes, feature_intensity_TICnorm_best_suffixes)
suffixes <- intersect(suffixes, feature_intensity_CLR_best_suffixes)
suffixes <- intersect(suffixes, feature_intensity_TICnorm_sum_suffixes)
suffixes <- intersect(suffixes, feature_count_suffixes)

# Step 2: Initialize a list to store results for each suffix
results_list <- list()
print(suffixes)
# Step 3: Process each suffix
for (suffix in suffixes) {

dt_process = copy(dt_tree_ids_redu)

    # Select relevant columns and rename to standard temporary names
    dt_temp <- dt_process[, .(
      FeatureID,
      UBERONBodyPartName,
      Cosine_temp = get(paste0("Cosine_", suffix)),
      Matching_Peaks_temp = get(paste0("Matching Peaks_", suffix)),
      inchi_key_first_block_count_temp = get(paste0("inchi_key_first_block_count_by_usi_", suffix)),
      feature_intensity_TICnorm_best_temp = get(paste0("feature_intensity_TICnorm_best_", suffix)),
      feature_intensity_CLR_best_temp = get(paste0("feature_intensity_CLR_best_", suffix)),
      feature_intensity_TICnorm_sum_temp = get(paste0("feature_intensity_TICnorm_sum_", suffix)),
      feature_count_temp = get(paste0("feature_count_", suffix))
    )]

dt_temp = dt_temp[Cosine_temp > 0]
    
    # Perform calculations on the renamed columns
    dt_temp <- dt_temp[, .(
        Cosine = max(Cosine_temp, na.rm = TRUE),
        feature_intensity_TICnorm_best = max(feature_intensity_TICnorm_best_temp, na.rm = TRUE),
        feature_intensity_CLR_best = max(feature_intensity_CLR_best_temp, na.rm = TRUE),
        feature_intensity_TICnorm_sum = max(feature_intensity_TICnorm_sum_temp, na.rm = TRUE),
        feature_count = max(feature_count_temp, na.rm = TRUE),
        Matching_Peaks = Matching_Peaks_temp[which.max(ifelse(is.na(Matching_Peaks_temp), -Inf, Matching_Peaks_temp))],
        inchi_key_first_block_count = {
            unique_values <- unique(inchi_key_first_block_count_temp[!is.na(inchi_key_first_block_count_temp)])
            if (length(unique_values) == 0) 0L else log10(length(unique(unlist(strsplit(as.character(unique_values), ",")))) + 0.1)
        },
        inchi_key_first_block_list = {
            unique_values <- unique(inchi_key_first_block_count_temp[!is.na(inchi_key_first_block_count_temp)])
            paste0(unique(unlist(strsplit(as.character(unique_values), ","))), collapse = ", ")
        }
    ), by = .(FeatureID, UBERONBodyPartName)]
    
    # Rename columns back to specific names with suffix
    setnames(dt_temp, c("Cosine", "Matching_Peaks", "inchi_key_first_block_count", "inchi_key_first_block_list", 
              "feature_intensity_TICnorm_best", "feature_intensity_CLR_best", "feature_intensity_TICnorm_sum", 
              "feature_count"),
             c(paste0("Cosine_", suffix), paste0("Matching_Peaks_", suffix), 
             paste0("inchi_key_first_block_count_by_usi_", suffix), 
               paste0("inchi_key_first_block_list_", suffix), 
               paste0("feature_intensity_TICnorm_best_", suffix), 
               paste0("feature_intensity_CLR_best_", suffix), 
               paste0("feature_intensity_TICnorm_sum_", suffix), 
               paste0("feature_count_", suffix)))

               print(suffix)

              print(nrow(dt_temp))    
    # Store result in the list
    results_list[[suffix]] <- dt_temp
}


print('Step 4')
# Step 4: Merge all results tables by FeatureID and UBERONBodyPartName
final_result <- Reduce(function(dt1, dt2) merge(dt1, dt2, by = c("FeatureID", "UBERONBodyPartName"), all = TRUE), results_list)


cols = c('FeatureID', 'UBERONBodyPartName',
         'tax_name', 'NCBIDivision', 'biotype', 
         'clade', 'class', 'cohort', 'family', 'genus', 'infraclass', 
         'infraorder', 'isolate', 'kingdom', 'order', 'parvorder', 
         'phylum', 'section', 'species', 'species group', 
         'species subgroup', 'strain', 'subclass', 'subcohort', 
         'subfamily', 'subgenus', 'subkingdom', 'suborder', 
         'subphylum', 'subsection', 'subspecies', 'subtribe', 
         'superclass', 'superfamily', 'superkingdom', 
         'superorder', 'tribe', 'varietas', 'NCBI')


# Add all columns starting with "sparql_"
sparql_cols <- grep("^sparql_", colnames(dt_tree_ids_redu), value = TRUE)
cols <- c(cols, sparql_cols)

print('step 5')

# print column named for both tables
print('final_result')
print(colnames(final_result))
print('\n\n\n\n\ndt_tree_ids_redu')
print(colnames(dt_tree_ids_redu))

# Step 5: Add back any other columns needed after merging
dt_tree_ids_redu <- merge(final_result, dt_tree_ids_redu[, ..cols], 
                           by = c("FeatureID", "UBERONBodyPartName"), all.y = TRUE, all.x = FALSE)

  print('done with this')
  #dt_databases = unique(dt_redu[Database != '' & !is.na(Database), c('uid_leaf', 'Database', 'NCBIDivision')])
  print(colnames(dt_tree_ids_redu))

  #fwrite(dt_tree_ids_redu[!is.na(uid_leaf) & uid_leaf != '', c('FeatureID', 'tax_name', 'detected', 'Cosine', 'Matching Peaks', 'UBERONBodyPartName', 'NCBIDivision')], paste0(c('masstResults_',  lib_id, '_', cid, '.tsv'), collapse = ''), sep = '\t')

  } else {
    dt_tree_ids_redu[, FeatureID := uid_leaf]
    dt_tree_ids_redu[, Cosine := NA]
    dt_tree_ids_redu[, detected := FALSE]
    dt_tree_ids_redu[, `Matching Peaks` := NA]
    dt_tree_ids_redu[, tax_name_class := NA]
    dt_tree_ids_redu[, tax_name_class_spec := NA]
    dt_tree_ids_redu[, inchi_key_first_block_count := NA]
    dt_tree_ids_redu[, inchi_key_first_block_list := NA]
    dt_tree_ids_redu[, feature_intensity_TICnorm_best := NA]
    dt_tree_ids_redu[, feature_intensity_CLR_best := NA]
    dt_tree_ids_redu[, feature_intensity_TICnorm_sum := NA]
    dt_tree_ids_redu[, feature_count := NA]
}

# Specify the fixed columns
fixed_cols <- c("FeatureID", "tax_name")

# Identify the stats columns
stats_cols <- setdiff(colnames(dt_tree_ids_redu), fixed_cols)

# Order the stats columns alphabetically
stats_cols_ordered <- stats_cols[order(stats_cols)]

# Reorder the data.table with the specified order
setcolorder(dt_tree_ids_redu, unique(c(fixed_cols, stats_cols_ordered)))


#set columns with substring inchi_key_first_block_count_by_usi_ 0 if they are NA
print('attemp to set 0 inchi_key_first_block_count_by_usi_')
inchi_key_cols = grep("inchi_key_first_block_count_by_usi_", colnames(dt_tree_ids_redu), value = TRUE)
dt_tree_ids_redu[, (inchi_key_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = inchi_key_cols]

print('attemp to set 0 sparql')
# sparql_cols = grep("sparql_", colnames(dt_tree_ids_redu), value = TRUE)
# dt_tree_ids_redu[, (sparql_cols) := lapply(.SD, function(x) fifelse(x==0, NA, x)), .SDcols = sparql_cols]
sparql_cols <- grep("sparql_", names(dt_tree_ids_redu), value = TRUE)
for (col in sparql_cols) {
  print(col)
  dt_tree_ids_redu[get(col) == 0, (col) := NA]
  print(paste(col, ' done'))
}



print('done')

dt_filtered = copy(dt_tree_ids_redu[!is.na(FeatureID) & FeatureID != '' & !grepl('mrcaott', FeatureID)])

print('filtered')

dt_output = unique(dt_filtered, by = 'FeatureID')

print('copied')

fwrite(dt_output, 'metadata_table_final.tsv', sep = '\t')





