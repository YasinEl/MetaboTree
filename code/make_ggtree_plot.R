# Load required libraries
library(optparse)
library(ggtree)
library(data.table)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
library(ape)
library(ggnewscale)

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
    } else if (any(grepl("Rosales", taxa_vector, ignore.case = TRUE))) {
      return("Rosales")
    } else if (any(grepl("Malvales", taxa_vector, ignore.case = TRUE))) {
      return("Malvales")
    } else if (any(grepl("Caryophyllales", taxa_vector, ignore.case = TRUE))) {
      return("Caryophyllales")
    } else if (any(grepl("Solanales", taxa_vector, ignore.case = TRUE))) {
      return("Solanales")
    } else if (any(grepl("Asparagales", taxa_vector, ignore.case = TRUE))) {
      return("Asparagales")
    } else if (any(grepl("Cytophagales", taxa_vector, ignore.case = TRUE))) {
      return("Cytophagales")
    } else if (any(grepl("Bryophyte", taxa_vector, ignore.case = TRUE))) {
      return("Bryophyte")
    } else if (any(grepl("Ranunculales", taxa_vector, ignore.case = TRUE))) {
      return("Ranunculales")     
    } else {
      return("Other")
    }
    # Fungi categories
  } else if (any(grepl("Xylariales", taxa_vector, ignore.case = TRUE))) {
    return("Xylariales")
  } else if (any(grepl("Helotiales", taxa_vector, ignore.case = TRUE))) {
    return("Helotiales")
  } else if (any(grepl("Hypocreales", taxa_vector, ignore.case = TRUE))) {
    return("Hypocreales")
  } else if (any(grepl("Polyporales", taxa_vector, ignore.case = TRUE))) {
    return("Polyporales")
    

    # Animal categories
  } else if (any(grepl("Animalia|Metazoa", taxa_vector, ignore.case = TRUE))) {
    if (any(grepl("Mammalia", taxa_vector, ignore.case = TRUE))) {
      return("Mammal")
    } else if (any(grepl("Aves", taxa_vector, ignore.case = TRUE))) {
      return("Bird")
    # } else if (any(grepl("Tunicata", taxa_vector, ignore.case = TRUE))) {
    #   return("Invertebrate")
    } else if (any(grepl("Reptilia|Serpentes", taxa_vector, ignore.case = TRUE))) {
      return("Reptile")
    } else if (any(grepl("Amphibia|Anura|Caudata|Gymnophiona", taxa_vector, ignore.case = TRUE))) {
      return("Amphibian")
    } else if (any(grepl("Pisces|Fish|Chondrichthyes|Osteichthyes|Actinopterygii|Sarcopterygii", taxa_vector, ignore.case = TRUE))) {
      return("Fish")
    } else if (any(grepl("Opisthokonta", taxa_vector, ignore.case = TRUE))) {
      return("Opisthokonta")


      
    } else {
      return("Other")
    }
    # Bacteria categories
  } else if (any(grepl("Bacteria", taxa_vector, ignore.case = TRUE))) {
    if (any(grepl("Pseudomonadota", taxa_vector, ignore.case = TRUE))) {
      return("Pseudomonadota")
    } else if (any(grepl("Actinomycetota", taxa_vector, ignore.case = TRUE))) {
      return("Actinomycetota")
    } else if (any(grepl("Bacillota", taxa_vector, ignore.case = TRUE))) {
      return("Bacillota")
    # } else if (any(grepl("Actinobacteria", taxa_vector, ignore.case = TRUE))) {
    #   return("Actinobacteria")
    } else if (any(grepl("Cyanobacteria", taxa_vector, ignore.case = TRUE))) {
      return("Cyanobacteria")
    } else if (any(grepl("Streptomyces", taxa_vector, ignore.case = TRUE))) {
      return("Streptomyces")

      
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
  make_option(c("-o", "--output_png"), type="character", default="tree.png", help="Output PNG file path", metavar="character")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

# Load tree
tree <- read.newick(args$input_tree)

# Load ReDU table
dt_redu <- fread(args$input_redu)

# Get Tree ID base table
dt_tree_ids = data.table(uid_leaf = tree$tip.label)

# Load data tables
dt_masst <- fread(args$input_masst)
if(nrow(dt_annotations) > 0){
  dt_masst[, match_col := USI2MASST_matchCol(unique(USI)), by =.(USI)]


  dt_masst[, max_cosine := max(Cosine), by = match_col]
  max_indices <- dt_masst[, .I[which.max(Cosine)], by = match_col]$V1
  dt_masst <- dt_masst[max_indices]

  df_merged <- merge(dt_masst, df_redu, by = "match_col", all = FALSE)
  print(paste("Merged MASST with REDU data with", nrow(df_merged), "rows."))

  df_merged[, n_detected := sum(!is.na(Cosine))]
  df_merged[, n_samples := .N]
  df_merged[, p_detected := 100 * n_detected / n_samples]



  dt_masst = dt_redu[dt_masst, on =.(match_col)]
  dt_masst = df_masst[uid_leaf %in% unique(dt_tree_ids$uid_leaf)]
  } else {
  dt_masst = data.table(uid_leaf = c(), p_detected = c())
}



# Prune tree
tree <- keep.tip(tree, intersect(unique(dt_redu$ID), tree$tip.label))







dt_ring_tissue = data.table(ID = tree$tip.label)


dt_redu[, bodypart := as.character(UBERONBodyPartName)]
dt_redu = dt_redu[, c('uid_leaf', 'bodypart')]

dt_redu[, n_bodypart := length(unique(ID)), by =.(bodypart)]

dt_redu[, n_bodyparts := length(unique(bodypart)), by =.(ID)]

dt_redu[n_bodyparts == 1 & n_bodypart >= 20, tissue_type := bodypart]
dt_redu[n_bodyparts == 1 & n_bodypart < 20, tissue_type := 'rare bodypart']
dt_redu[n_bodyparts > 2, tissue_type := 'multiple']
dt_redu[n_bodyparts == 2, tissue_type := {
  if('missing value' %in% bodypart){
    bp = unique(bodypart)[unique(bodypart) != 'missing value']
    rarity = unique(n_bodypart)[unique(bodypart) != 'missing value']
    if (rarity < 20){
      'rare bodypart'
    } else {
      bp
    }
  } else {
    'multiple'
  }
  
}, by =.(ID)]


dt_bodypart = unique(dt_redu[, c('ID', 'tissue_type')])
dt_bodypart[tissue_type != 'gallbladder', tissue_type := 'Other']


dt_ring_tissue = dt_bodypart[dt_ring_tissue, on =.(ID)]

if(nrow(dt_annotations) > 0){
dt_annotations = dt_bodypart[dt_annotations, on =.(ID)]
}


print(dt_annotations)

#source
dt_redu = fread(args$input_redu)
dt_redu = dt_redu[grepl('|', NCBITaxonomy, fixed = TRUE)]
dt_redu[, ncbiid := strsplit(unique(NCBITaxonomy), '|', fixed = TRUE)[[1]][1], by =.(NCBITaxonomy)]
dt_redu[, ID := as.character(ncbiid)]

dtDataSource = dt_redu[, .(DataSource = ifelse(length(unique(DataSource)) > 1, 'multiple', unique(DataSource))) , by =.(ID)]



dt_sparql = fread(args$input_sparql)
if(nrow(dt_sparql) > 0) {
dt_sparql[, ID := as.character(ID)]
dt_sparql[, node := ID]
dt_sparql = dt_sparql[, c('node', 'present')]
dt_sparql[, node := as.integer(node)]
dt_sparql[, x_val := 1]
}





dt_ncbi = fread(args$input_ncbi)
if(nrow(dt_ncbi) > 0){
dt_ncbi[, ID := as.character(ID)]
dt_ncbi[, node := ID]
dt_ncbi = dt_ncbi[, c('node', 'present')]
dt_ncbi[, node := as.integer(node)]
dt_ncbi[, x_val := 1]
}






p_t = 
ggtree(tree, layout='circular', size=0.15, open.angle=5) + 
  geom_fruit(data=dt_ring_kingdom, geom=geom_tile,
             mapping=aes(y=ID,  fill=kingdom),
             pwidth = 1)  +
  scale_fill_manual(values=c(
    "Bacteria" = "#17becf",
    "Animal" = "#ff7f0e",
    "Fungi" = "#9467bd",
    "Plant" = "#2ca02c",
    "Other" = "#7f7f7f",
    "Archaea" = "#1f77b4"
  )) +
  new_scale_fill() +
  geom_fruit(data=dt_ring_superclass, geom=geom_tile,
             mapping=aes(y=ID,  fill=superclass),
             pwidth = 1) +
  scale_fill_manual(values=
    category_colors <- c(     
      "Proteobacteria" = "#8A2BE2",      # Blue Violet: vibrant for diverse bacteria group
      "Pseudomonadota" = "#DEB887",          # Burly Wood: earthy tone for robust bacteria
      "Actinomycetota" = "#A52A2A",      # Brown: earthy, representing soil habitats
      "Fish" = "#4682B4",                # Steel Blue: aquatic color, typical for fish
      "Mammal" = "#A52A2D",              # Sienna: robust and earthy
      "Bird" = "#F0E68C",                # Khaki: airy and light, like feathers
      "Amphibian" = "#556B2F",           # Pink: light and varied for miscellaneous animals
      "Invertebrate" = "#FFD700",        # Gold: diverse and fascinating
      "Echinoderm" = "#FF6347",          # Tomato: vibrant, oceanic organisms
      "Crustacean" = "#FA8072",          # Salmon: related to their often reddish hue
      "Insect" = "#ADFF2F",              # Green Yellow: vibrant, representing diversity
      "Arachnid" = "#808000",            # Olive: dark and mysterious
      "Nematode" = "#FF4500",            # Orange Red: vivid, alerting to their presence
      "Mollusk" = "#B0E0E6",             # Powder Blue: marine and fresh, like many mollusks
      "Annelid" = "#9932CC",             # Dark Orchid: earthy and underexplored
      "Cnidarian" = "#1E90FF",           # Dodger Blue: oceanic and vibrant
      "Sponge" = "#9ACD32",              # Yellow Green: simple, foundational marine life
      "Bacillota" = "#DAA520",               # Sea Green: generic vibrant plant color
      "Sapindales" = "#FF4500",          # Orange Red: vibrant, lively plant group
      "Malpighiales" = "#3CB371",        # Medium Sea Green: lush, leafy
      "Fabales" = "#20B2AA",             # Light Sea Green: fresh, leguminous
      "Rosales" = "#DB7093",             # Pale Violet Red: floral and delicate
      "Actinobacteria" = "#6A5ACD",           # Gold: bright and impactful
      "Lamiales" = "#9ACD32",            # Yellow Green: common, everyday flora
      "Gentianales" = "#40E0D0",         # Turquoise: medicinal and tropical
      "Caryophyllales" = "#EE82EE",      # Violet: vibrant and varied
      "Poales" = "#48D1CC",              # Medium Turquoise: grassy and widespread
      "Other" = "#C0C0C0",               # Silver: neutral for unspecified categories
      "Archaea" = "#87CEFA",              # Light Sky Blue: ancient and fundamental
      "Streptomyces" = "#A52A2D",
      "Actinomycetota" = "#48D1CC",
      "Eurotiomycetidae" = "#1E90FF",
      "Cyanobacteria" = "#C0C0C0",
      "Streptomyces" = "#DAA520",
      "Asparagales" = "#4682B4",
      "Malvales" = "#FFD700",
#Plant          
      "Solanales" = "#FFD700", 
      "Apiales" = "#FA8072",
      "Ranunculales" = "#40E0D0",
#Fungi
      "Xylariales" = "#FFD700",
      "Helotiales" = "#FA8072",
#Animal
      "Reptile" = "#3CB371",
      "Opisthokonta" = "#9ACD32"
    )

  ) +
  # new_scale_fill() +
  # geom_fruit(data=dt_ring_tissue, geom=geom_tile,
  #            mapping=aes(y=ID, fill=tissue_type),
  #            pwidth = 1) +
  # scale_fill_manual(values = c(
  #   "multiple" = "#6c71c4",        # A cool lavender color, suggests variety
  #   "rare bodypart" = "#cb4b16",   # A bold red-orange, indicating rarity
  #   "missing value" = "#93a1a1",   # Grey, often used for missing or unavailable data
  #   "feces" = "#586e75",           # Dark slate, earthy and organic
  #   "gallbladder" = "#b58900",     # A deep yellow, reflective of bile
  #   "digestive tract" = "#dc322f", # Red, associated with the internal organ color
  #   "leaf" = "#859900",            # A vibrant green, typical for leaves
  #   "root" = "#b15928"             # Brownish-orange, reminiscent of soil and roots
  # )) +
  # new_scale_fill() +
  # geom_fruit(data=dtDataSource, geom=geom_tile,
  #            mapping=aes(y=ID, fill=DataSource),
  #            pwidth = 1) +
  # scale_fill_manual(values = c(
  #   "multiple" = "#6c71c4",        # A cool lavender color, suggests variety
  #   "GNPS" = "#cb4b16",   # A bold red-orange, indicating rarity
  #   "Workbench" = "#FFD700",   # Grey, often used for missing or unavailable data
  #   "MetaboLights" = "#586e75"           # Brownish-orange, reminiscent of soil and roots
  # )) +
  new_scale_fill() 


if(nrow(dt_sparql) > 0){

  p_t = p_t + geom_fruit(data=dt_sparql, geom_point,
  mapping=aes(y=node), size=1, alpha=0.5, color = 'red') +
  new_scale_fill()

}

if(nrow(dt_ncbi) > 0){

  p_t = p_t + geom_fruit(data=dt_ncbi, geom_point,
  mapping=aes(y=node), size=1, alpha=0.5, color = 'green') +
  new_scale_fill()

}

if(nrow(dt_annotations) > 0){
  p_t = 
  p_t + geom_fruit(data=dt_annotations, geom=geom_bar,
            mapping=aes(y=ID, x=p_detected, fill = tissue_type),
            width = 1,
            pwidth=0.5, 
            orientation="y", 
            stat="identity"
            ) +
  scale_fill_manual(values = c(
    "multiple" = "#6c71c4",        # A cool lavender color, suggests variety
    "rare bodypart" = "#cb4b16",   # A bold red-orange, indicating rarity
    "missing value" = "#93a1a1",   # Grey, often used for missing or unavailable data
    "feces" = "#586e75",           # Dark slate, earthy and organic
    "gallbladder" = "#b58900",     # A deep yellow, reflective of bile
    "digestive tract" = "#dc322f", # Red, associated with the internal organ color
    "leaf" = "#859900",            # A vibrant green, typical for leaves
    "root" = "#b15928",             # Brownish-orange, reminiscent of soil and roots
    "Other" = "#87CEFA"            
  )) #+
  # guides(
  #   fill=guide_legend(ncol=1, title="Legend"),
  #   fill.1=guide_legend(ncol=1, title="Superclass"),
  #   fill.2=guide_legend(ncol=1, title="Tissue Type"),
  #   fill.3=guide_legend(ncol=2, title="Data Source"),
  #   fill.4=guide_legend(ncol=1, title="Annotations")
  # )


}


# p_t = p_t %<+% geom_nodepoint(dt_sparql, aes(subset=node))



ggsave(args$output_png, plot = p_t, width = 10, height = 8, dpi = 300)