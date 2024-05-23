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


classify_taxa <- function(taxa_vector) {
  # Check each category with their respective patterns
  if (any(grepl("Plantae", taxa_vector, ignore.case = TRUE))) {
    return("Plant")
  } else if (any(grepl("Fungi", taxa_vector, ignore.case = TRUE))) {
    return("Fungi")
  } else if (any(grepl("Animalia|Metazoa", taxa_vector, ignore.case = TRUE))) {
    return("Animal")
  } else if (any(grepl("Bacteria", taxa_vector, ignore.case = TRUE))) {
    return("Bacteria")
  } else if (any(grepl("Archaea", taxa_vector, ignore.case = TRUE))) {
    return("Archaea")  # Added condition for Archaea
  } else {
    return("Other")
  }
}
classify_taxa_detailed <- function(taxa_vector) {
  # Plant categories
  if (any(grepl("Plantae", taxa_vector, ignore.case = TRUE))) {
    if (any(grepl("Lamiales", taxa_vector, ignore.case = TRUE))) {
      return("Lamiales")
    } else if (any(grepl("Asterales", taxa_vector, ignore.case = TRUE))) {
      return("Asterales")
    } else if (any(grepl("Fabales", taxa_vector, ignore.case = TRUE))) {
      return("Fabales")
    } else if (any(grepl("Malpighiales", taxa_vector, ignore.case = TRUE))) {
      return("Malpighiales")
    } else if (any(grepl("Rosales", taxa_vector, ignore.case = TRUE))) {
      return("Rosales")
    } else if (any(grepl("Gentianales", taxa_vector, ignore.case = TRUE))) {
      return("Gentianales")
    } else if (any(grepl("Sapindales", taxa_vector, ignore.case = TRUE))) {
      return("Sapindales")
    } else if (any(grepl("Poales", taxa_vector, ignore.case = TRUE))) {
      return("Poales")
    } else if (any(grepl("Caryophyllales", taxa_vector, ignore.case = TRUE))) {
      return("Caryophyllales")
    } else if (any(grepl("Solanales", taxa_vector, ignore.case = TRUE))) {
      return("Solanales")
    } else {
      return("Other")
    }
    # Fungi categories
  } else if (any(grepl("Fungi", taxa_vector, ignore.case = TRUE))) {
    return("Fungi")
    # Animal categories
  } else if (any(grepl("Animalia|Metazoa", taxa_vector, ignore.case = TRUE))) {
    if (any(grepl("Mammalia", taxa_vector, ignore.case = TRUE))) {
      return("Mammal")
    } else if (any(grepl("Aves", taxa_vector, ignore.case = TRUE))) {
      return("Bird")
    } else if (any(grepl("Tunicata", taxa_vector, ignore.case = TRUE))) {
      return("Invertebrate")
    } else if (any(grepl("Reptilia", taxa_vector, ignore.case = TRUE))) {
      return("Reptile")
    } else if (any(grepl("Amphibia|Anura|Caudata|Gymnophiona", taxa_vector, ignore.case = TRUE))) {
      return("Amphibian")
    } else if (any(grepl("Pisces|Fish|Chondrichthyes|Osteichthyes|Actinopterygii|Sarcopterygii", taxa_vector, ignore.case = TRUE))) {
      return("Fish")
    } else if (any(grepl("Echinodermata|Echinacea|Echinidea|Echinoidea|Euechinoidea", taxa_vector, ignore.case = TRUE))) {
      return("Echinoderm")
    } else if (any(grepl("Nematoda|Anisakidae|Ascaridoidea|Ascaridomorpha", taxa_vector, ignore.case = TRUE))) {
      return("Nematode")
    } else if (any(grepl("Cnidaria|Alcyonium|Plexaura", taxa_vector, ignore.case = TRUE))) {
      return("Cnidarian")
    } else if (any(grepl("Annelida|Parechinidae|Polychaeta|Terebellida|Terebellidae|Terebelliformia", taxa_vector, ignore.case = TRUE))) {
      return("Annelid")
    } else if (any(grepl("Insecta", taxa_vector, ignore.case = TRUE))) {
      return("Insect")
    } else if (any(grepl("Crustacea", taxa_vector, ignore.case = TRUE))) {
      return("Crustacean")
    } else if (any(grepl("Mollusca", taxa_vector, ignore.case = TRUE))) {
      return("Mollusk")
    } else if (any(grepl("Arachnida", taxa_vector, ignore.case = TRUE))) {
      return("Arachnid")
    } else if (any(grepl("Chondrichthyes", taxa_vector, ignore.case = TRUE))) {
      return("Cartilaginous Fish")
    } else if (any(grepl("Actinopterygii", taxa_vector, ignore.case = TRUE))) {
      return("Ray-finned Fish")
    } else if (any(grepl("Serpentes", taxa_vector, ignore.case = TRUE))) {
      return("Snake")
    } else if (any(grepl("Porifera", taxa_vector, ignore.case = TRUE))) {
      return("Sponge")
    } else if (any(grepl("Cephalopoda", taxa_vector, ignore.case = TRUE))) {
      return("Cephalopod")
    } else if (any(grepl("Gastropoda", taxa_vector, ignore.case = TRUE))) {
      return("Gastropod")
    } else {
      return("Other")
    }
    # Bacteria categories
  } else if (any(grepl("Bacteria", taxa_vector, ignore.case = TRUE))) {
    if (any(grepl("Proteobacteria", taxa_vector, ignore.case = TRUE))) {
      return("Proteobacteria")
    } else if (any(grepl("Firmicutes", taxa_vector, ignore.case = TRUE))) {
      return("Firmicutes")
    } else if (any(grepl("Actinobacteria", taxa_vector, ignore.case = TRUE))) {
      return("Actinobacteria")
    } else if (any(grepl("Cyanobacteria", taxa_vector, ignore.case = TRUE))) {
      return("Cyanobacteria")
    } else {
      return("Other")
    }
    # Archaea categories
  } else if (any(grepl("Archaea", taxa_vector, ignore.case = TRUE))) {
    return("Archaea")
  } else {
    return("Other")
  }
}



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
  make_option(c("-u", "--usi"), type="character", default=NULL, help="USI", metavar="character"),
  make_option(c("-p", "--cid"), type="character", default=NULL, help="CID", metavar="character")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

lib_id = strsplit(args$usi, ':')[[1]][5]
cid = args$cid

# Load ReDU table
dt_redu <- fread(args$input_redu)
#dt_redu[, n_samples := sum(!is.na(match_col)), by =.(uid_leaf)]

dt_redu_match = dt_redu[!is.na(uid_leaf), .(NCBIDivision = NCBIDivision[1],UBERONBodyPartName = paste0(sort(unique(UBERONBodyPartName)), collapse= ', '), tax_name = tax_name[1], Database = unique(Database)[1]), by =.(uid_leaf)]

# Load tree
tree <- read.newick(args$input_tree)

# Get Tree ID base table
dt_tree_ids = data.table(uid_leaf = c(tree$tip.label, tree$node.label))

print(head(dt_redu_match))
print(colnames(dt_redu_match))

dt_tree_ids_redu <- merge(dt_tree_ids, dt_redu_match, by = "uid_leaf", all.x = TRUE)

fwrite(dt_tree_ids_redu, 'check_dt_tree_ids_redu.csv')

# Load data tables
dt_masst <- fread(args$input_masst)

if(nrow(dt_masst) > 0){
  dt_masst[, match_col := USI2MASST_matchCol(unique(USI)), by =.(USI)]


  dt_masst[, max_cosine := max(Cosine), by = match_col]
  max_indices <- dt_masst[, .I[which.max(Cosine)], by = match_col]$V1
  dt_masst <- dt_masst[max_indices]
  dt_masst <- unique(dt_masst[, c('Cosine', 'Matching Peaks', 'match_col')])

  dt_masst <- dt_masst[, .(Cosine = max(Cosine), `Matching Peaks` = max(`Matching Peaks`[which.max(Cosine)])), by = .(match_col)]

  print('dt_tree_ids_redu')
  print(head(dt_tree_ids_redu))
  print(colnames(dt_tree_ids_redu))

  print('dt_masst')
  print(head(dt_masst))
  print(colnames(dt_masst))

  dt_redu_match = unique(dt_redu[, c('match_col', 'uid_leaf')])


  fwrite(dt_masst, 'check_masst.csv')

  dt_redu_match <- merge(dt_masst, dt_redu_match, by = "match_col", all = TRUE)

  dt_redu_match = unique(dt_redu_match)

  fwrite(dt_redu_match, 'dt_redu_match.csv')

  dt_tree_ids_redu <- merge(dt_tree_ids_redu, dt_redu_match, by = "uid_leaf", all = TRUE)

  dt_tree_ids_redu[, detected := FALSE]
  dt_tree_ids_redu[!is.na(Cosine), detected := TRUE]
  dt_tree_ids_redu[, FeatureID := uid_leaf]

  dt_tree_ids_redu = dt_tree_ids_redu[, .(Cosine = max(Cosine, na.rm = TRUE), `Matching Peaks` = `Matching Peaks`[which.max(`Matching Peaks`)], detected = any(detected)), by =.(FeatureID, tax_name, UBERONBodyPartName, NCBIDivision, Database)]

  #dt_databases = unique(dt_redu[Database != '' & !is.na(Database), c('uid_leaf', 'Database', 'NCBIDivision')])


  dt_tree_ids_redu[, tax_name_class := classify_taxa(unique(tax_name)), by =.(tax_name)]
  dt_tree_ids_redu[, tax_name_class_spec := classify_taxa_detailed(unique(tax_name)), by =.(tax_name)]


  #fwrite(dt_tree_ids_redu[!is.na(uid_leaf) & uid_leaf != '', c('FeatureID', 'tax_name', 'detected', 'Cosine', 'Matching Peaks', 'UBERONBodyPartName', 'NCBIDivision')], paste0(c('masstResults_',  lib_id, '_', cid, '.tsv'), collapse = ''), sep = '\t')

  } else {
    dt_tree_ids_redu[, FeatureID := uid_leaf]
    dt_tree_ids_redu[, Cosine := NA]
    dt_tree_ids_redu[, detected := FALSE]
    dt_tree_ids_redu[, `Matching Peaks` := NA]
    dt_tree_ids_redu[, tax_name_class := NA]
    dt_tree_ids_redu[, tax_name_class_spec := NA]
}

fwrite(dt_tree_ids_redu[!is.na(FeatureID) & FeatureID != '', c('FeatureID', 'tax_name', 'detected', 'Cosine', 'Matching Peaks', 'UBERONBodyPartName', 'NCBIDivision', 'Database', 'tax_name_class', 'tax_name_class_spec')], paste0(c('treeAnnotation_',  lib_id, '_', cid, '.tsv'), collapse = ''), sep = '\t')





