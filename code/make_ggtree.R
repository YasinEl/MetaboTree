# Load required libraries
library(optparse)
library(ggtree)
library(data.table)
library(ggplot2)



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



# Setup command-line options
option_list <- list(
  make_option(c("-i", "--input_tree"), type="character", default=NULL, help="Input Newick file path", metavar="character"),
  make_option(c("-l", "--input_lib"), type="character", default=NULL, help="Input library file path", metavar="character"),
  make_option(c("-lin", "--input_lin"), type="character", default=NULL, help="Input linage file path", metavar="character"),
  make_option(c("-r", "--input_redu"), type="character", default=NULL, help="Input REDU file path", metavar="character"),
  make_option(c("-o", "--output_png"), type="character", default="tree.png", help="Output PNG file path", metavar="character")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

# Load tree
tree <- read.newick(args$input_tree)

# Load data tables
dt_library <- fread(args$input_lib)
dt_redu <- fread(args$input_redu)
dt_redu <- dt_redu[grepl('|', NCBITaxonomy, fixed = TRUE)]
dt_redu[, ncbiid := strsplit(unique(NCBITaxonomy), '|', fixed = TRUE)[[1]][1], by = .(NCBITaxonomy)]
dt_redu[, ID := as.character(ncbiid)]

# Prune tree
tree <- keep.tip(tree, intersect(unique(dt_redu$ID), tree$tip.label))



dt_ring_kingdom = data.table(ID = tree$tip.label)

dt_redu_linages = fread(args$input_lin)
dt_redu_linages[, kingdom := classify_taxa(name), by =.(ncbi)]
dt_redu_linages[, ID := as.character(ncbi)]
dt_redu_linages = unique(dt_redu_linages[, c('ID', 'kingdom')])
dt_redu_linages = dt_redu_linages[ID %in% dt_ring_kingdom$ID]

dt_ring_kingdom = dt_redu_linages[dt_ring_kingdom, on =.(ID)]




dt_ring_superclass = data.table(ID = tree$tip.label)

dt_redu_linages = fread(args$input_lin)
dt_redu_linages[, superclass := classify_taxa_detailed(name), by =.(ncbi)]
dt_redu_linages[, ID := as.character(ncbi)]
dt_redu_linages = unique(dt_redu_linages[, c('ID', 'superclass')])
dt_redu_linages = dt_redu_linages[ID %in% dt_ring_superclass$ID]

dt_ring_superclass = dt_redu_linages[dt_ring_superclass, on =.(ID)]


dt_ring_tissue = data.table(ID = tree$tip.label)


dt_redu[, bodypart := as.character(UBERONBodyPartName)]
dt_redu = dt_redu[, c('ID', 'bodypart')]

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


dt_ring_tissue = dt_bodypart[dt_ring_tissue, on =.(ID)]



#source
dt_redu = fread(args$input_redu)
dt_redu = dt_redu[grepl('|', NCBITaxonomy, fixed = TRUE)]
dt_redu[, ncbiid := strsplit(unique(NCBITaxonomy), '|', fixed = TRUE)[[1]][1], by =.(NCBITaxonomy)]
dt_redu[, ID := as.character(ncbiid)]

dtDataSource = dt_redu[, .(DataSource = ifelse(length(unique(DataSource)) > 1, 'multiple', unique(DataSource))) , by =.(ID)]



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
      "Other Bacteria" = "#C0C0C0",      # Light Steel Blue: neutral for miscellaneous bacteria
      "Proteobacteria" = "#8A2BE2",      # Blue Violet: vibrant for diverse bacteria group
      "Firmicutes" = "#DEB887",          # Burly Wood: earthy tone for robust bacteria
      "Actinobacteria" = "#A52A2A",      # Brown: earthy, representing soil habitats
      "Fish" = "#4682B4",                # Steel Blue: aquatic color, typical for fish
      "Mammal" = "#A52A2D",              # Sienna: robust and earthy
      "Bird" = "#F0E68C",                # Khaki: airy and light, like feathers
      "Amphibian" = "#556B2F",           # Dark Olive Green: wetland habitats
      "Other Animal" = "#C0C0C0",        # Pink: light and varied for miscellaneous animals
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
      "Fungi" = "#DAA520",               # Golden Rod: earthy and fungal
      "Other Plant" = "#C0C0C0",         # Sea Green: generic vibrant plant color
      "Sapindales" = "#FF4500",          # Orange Red: vibrant, lively plant group
      "Malpighiales" = "#3CB371",        # Medium Sea Green: lush, leafy
      "Fabales" = "#20B2AA",             # Light Sea Green: fresh, leguminous
      "Rosales" = "#DB7093",             # Pale Violet Red: floral and delicate
      "Asterales" = "#6A5ACD",           # Slate Blue: diverse and widespread
      "Solanales" = "#FFD700",           # Gold: bright and impactful
      "Lamiales" = "#9ACD32",            # Yellow Green: common, everyday flora
      "Gentianales" = "#40E0D0",         # Turquoise: medicinal and tropical
      "Caryophyllales" = "#EE82EE",      # Violet: vibrant and varied
      "Poales" = "#48D1CC",              # Medium Turquoise: grassy and widespread
      "Other" = "#C0C0C0",               # Silver: neutral for unspecified categories
      "Archaea" = "#87CEFA"              # Light Sky Blue: ancient and fundamental
    )
    
  ) +
  new_scale_fill() +
  geom_fruit(data=dt_ring_tissue, geom=geom_tile,
             mapping=aes(y=ID, fill=tissue_type),
             pwidth = 1) +
  scale_fill_manual(values = c(
    "multiple" = "#6c71c4",        # A cool lavender color, suggests variety
    "rare bodypart" = "#cb4b16",   # A bold red-orange, indicating rarity
    "missing value" = "#93a1a1",   # Grey, often used for missing or unavailable data
    "feces" = "#586e75",           # Dark slate, earthy and organic
    "gallbladder" = "#b58900",     # A deep yellow, reflective of bile
    "digestive tract" = "#dc322f", # Red, associated with the internal organ color
    "leaf" = "#859900",            # A vibrant green, typical for leaves
    "root" = "#b15928"             # Brownish-orange, reminiscent of soil and roots
  )) +
  new_scale_fill() +
  geom_fruit(data=dtDataSource, geom=geom_tile,
             mapping=aes(y=ID, fill=DataSource),
             pwidth = 1) +
  scale_fill_manual(values = c(
    "multiple" = "#6c71c4",        # A cool lavender color, suggests variety
    "GNPS" = "#cb4b16",   # A bold red-orange, indicating rarity
    "Workbench" = "#FFD700",   # Grey, often used for missing or unavailable data
    "MetaboLights" = "#586e75"           # Brownish-orange, reminiscent of soil and roots
  )) +
  new_scale_fill() +
  geom_fruit(data = dt_library, geom=geom_tile,
        mapping=aes(y=ID, x = molecule, fill=fill),
        pwidth = 1)


ggsave(args$output_png, plot = p_t, width = 10, height = 8, dpi = 300)