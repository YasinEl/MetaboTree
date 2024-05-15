# Load required libraries
library(optparse)
library(ggtree)
library(data.table)
library(ggplot2)
library(ggtreeExtra)
library(treeio)
library(ape)
library(ggnewscale)
library(stringr)
#library(phytools)


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
  make_option(c("-p", "--cid"), type="character", default=NULL, help="CID", metavar="character"),
  make_option(c("-o", "--output_png"), type="character", default="tree.png", help="Output PNG file path", metavar="character")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

lib_id = strsplit(args$usi, ':')[[1]][5]
cid = args$cid

# Load ReDU table
dt_redu <- fread(args$input_redu)
dt_redu[, n_samples := sum(!is.na(match_col)), by =.(uid_leaf)]

dt_ncbi_kingdom = unique(dt_redu[, c('uid_leaf', 'NCBIDivision')])


# Load tree
tree <- read.newick(args$input_tree)



# Get Tree ID base table
dt_tree_ids = data.table(uid_leaf = c(tree$tip.label, tree$node.label))


# Load data tables
dt_masst <- fread(args$input_masst)


if(nrow(dt_masst) > 0){
  dt_masst[, match_col := USI2MASST_matchCol(unique(USI)), by =.(USI)]


  dt_masst[, max_cosine := max(Cosine), by = match_col]
  max_indices <- dt_masst[, .I[which.max(Cosine)], by = match_col]$V1
  dt_masst <- dt_masst[max_indices]

  

  dt_masst <- merge(dt_masst, dt_redu, by = "match_col", all = FALSE)

  dt_databases = unique(dt_redu[Database != '' & !is.na(Database), c('uid_leaf', 'Database', 'NCBIDivision')])
  #dt_masst = rbindlist(list(dt_masst, dt_redu_extend), fill = TRUE, use.names = TRUE)

  print(paste("Merged MASST with REDU data with", nrow(dt_masst), "rows."))

fwrite(dt_databases, 'dt_databases.csv')

  dt_masst = dt_masst[uid_leaf %in% unique(dt_tree_ids$uid_leaf)]

fwrite(dt_masst, 'bdf_merged_check1.csv')

  dt_masst = dt_masst[, .(n_detected = sum(!is.na(match_col)), p_detected = sum(!is.na(match_col))/unique(n_samples)*100), by =.(uid_leaf, UBERONBodyPartName, NCBIDivision, Database)]

  fwrite(dt_masst, 'cdf_merged_processed.csv')

  } else {
  dt_masst = data.table(uid_leaf = c(), p_detected = c())
}





dt_masst[, p_detected := as.numeric(p_detected)]

p_t = ggtree(tree, layout='circular', size=0.15, open.angle=5) + 
geom_point2(aes(subset = (node == 'ott84004')), size = 10, color = "red") +
ggtitle(paste0(c('GNPS:', lib_id, 'Pubchem:', cid), collapse = ' ')) +
  geom_fruit(data=dt_ncbi_kingdom, geom=geom_tile,
             mapping=aes(y=uid_leaf,  fill=NCBIDivision),
             pwidth = 1)  +
  scale_fill_manual(values=c(
    "Bacteria" = "#1f77b4",          # Blue
    "Mammals" = "#ff7f0e",           # Orange
    "Primates" = "#d62728",          # Red
    "Invertebrates" = "#9467bd",     # Purple
    "Plants and Fungi" = "#2ca02c",  # Green
    "Rodents" = "#bcbd22",           # Olive
    "Vertebrates" = "#17becf",       # Cyan
    "missing value" = "#7f7f7f"      # Gray
  )) +
  new_scale_fill() 

if(nrow(dt_databases) > 0){

  p_t = p_t + geom_fruit(data=dt_databases, geom_point,
  mapping=aes(y=uid_leaf, color = Database), size=1, alpha=0.5) +
  new_scale_color()

}

if(nrow(dt_masst) > 0){
  p_t = 
  p_t + geom_fruit(data=dt_masst, geom=geom_point,
            mapping=aes(y=uid_leaf, color = p_detected)) +
    scale_color_continuous(name = "p_detected") #+
}


# if(nrow(dt_masst) > 0){
#   p_t = 
#   p_t + geom_fruit(data=dt_masst, geom=geom_bar,
#             mapping=aes(y=uid_leaf, x=p_detected, fill = UBERONBodyPartName),
#             width = 1,
#             pwidth=0.5, 
#             orientation="y", 
#             stat="identity"
#             ) #+
#   # scale_fill_manual(values = c(
#   #   "multiple" = "#6c71c4",        # A cool lavender color, suggests variety
#   #   "rare bodypart" = "#cb4b16",   # A bold red-orange, indicating rarity
#   #   "missing value" = "#93a1a1",   # Grey, often used for missing or unavailable data
#   #   "feces" = "#586e75",           # Dark slate, earthy and organic
#   #   "gallbladder" = "#b58900",     # A deep yellow, reflective of bile
#   #   "digestive tract" = "#dc322f", # Red, associated with the internal organ color
#   #   "leaf" = "#859900",            # A vibrant green, typical for leaves
#   #   "root" = "#b15928",             # Brownish-orange, reminiscent of soil and roots
#   #   "Other" = "#87CEFA"            
#   # )) #+
#   # guides(
#   #   fill=guide_legend(ncol=1, title="Legend"),
#   #   fill.1=guide_legend(ncol=1, title="Superclass"),
#   #   fill.2=guide_legend(ncol=1, title="Tissue Type"),
#   #   fill.3=guide_legend(ncol=2, title="Data Source"),
#   #   fill.4=guide_legend(ncol=1, title="Annotations")
#   # )


# }


# p_t = p_t %<+% geom_nodepoint(dt_sparql, aes(subset=node))



ggsave(paste0(c(lib_id, '_', cid, '.png'), collapse = ''), plot = p_t, width = 10, height = 8, dpi = 300)