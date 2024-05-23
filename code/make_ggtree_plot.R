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

save_plot_with_svg <- function(plot, svg_grob, filename) {
  # Open a graphics device
  png(filename, width = 2400, height = 1800, res = 300)
  
  # Print the ggplot
  print(plot)
  
  # Create a viewport for the subplot in the upper right corner
  vp <- viewport(width = 0.5, height = 0.5, x = 1, y = 0.9)
  
  # Overlay the SVG subplot
  pushViewport(vp)
  grid.draw(svg_grob)
  popViewport()
  
  # Close the graphics device
  dev.off()
}


# Setup command-line options
option_list <- list(
  make_option(c("-i", "--input_tree"), type="character", default=NULL, help="Input Newick file path", metavar="character"),
  make_option(c("-l", "--input_masst"), type="character", default=NULL, help="Input library file path", metavar="character"),
  make_option(c("-r", "--input_redu"), type="character", default=NULL, help="Input REDU file path", metavar="character"),
  make_option(c("-m", "--mol_plot"), type="character", default=NULL, help="Molecule png path", metavar="character"),
  make_option(c("-n", "--mol_name"), type="character", default=NULL, help="Molecule name tsv path", metavar="character"),
  make_option(c("-u", "--usi"), type="character", default=NULL, help="USI", metavar="character"),
  make_option(c("-p", "--cid"), type="character", default=NULL, help="CID", metavar="character"),
  make_option(c("-o", "--output_png"), type="character", default="tree.png", help="Output PNG file path", metavar="character")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

lib_id = strsplit(args$usi, ':')[[1]][5]
cid = args$cid


#load moelcule name
name_df <- read_tsv(args$mol_name, col_names = c("Name"))
molecule_name <- as.character(name_df$Name)



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

dt_databases = data.table(uid_leaf = c(), Database = c(), NCBIDivision = c())

if(nrow(dt_masst) > 0){
  dt_masst[, match_col := USI2MASST_matchCol(unique(USI)), by =.(USI)]


  dt_masst[, max_cosine := max(Cosine), by = match_col]
  max_indices <- dt_masst[, .I[which.max(Cosine)], by = match_col]$V1
  dt_masst <- dt_masst[max_indices]

  

  dt_masst <- merge(dt_masst, dt_redu, by = "match_col", all = FALSE)

  dt_databases = unique(dt_redu[Database != '' & !is.na(Database), c('uid_leaf', 'Database', 'NCBIDivision')])

  print(paste("Merged MASST with REDU data with", nrow(dt_masst), "rows."))

  dt_masst = dt_masst[uid_leaf %in% unique(dt_tree_ids$uid_leaf)]
  


  dt_masst_exp = dt_masst[, .(Cosine = max(Cosine )), by =.(uid_leaf, tax_name)]

  dt_masst_exp [, FeatureID := uid_leaf]

  fwrite(dt_masst_exp[!is.na(uid_leaf) & uid_leaf != '', c('FeatureID', 'tax_name', 'Cosine')], paste0(c('masstResults_',  lib_id, '_', cid, '.tsv'), collapse = ''), sep = '\t')

  dt_masst = dt_masst[!is.na(Cosine)]

  #fwrite(dt_masst[, c('uid_leaf', 'Database', 'tax_name', 'NCBI', 'NCBIDivision', 'UBERONBodyPartName', 'Delta Mass', 'Cosine', 'Matching Peaks', 'USI')], paste0(c('masstResults_',  lib_id, '_', cid, '.tsv'), collapse = ''), sep = '\t')


  dt_masst = dt_masst[, .(n_detected = sum(!is.na(match_col)), p_detected = sum(!is.na(match_col))/unique(n_samples)*100, Cosine = max(Cosine)), by =.(uid_leaf)]
  } else {
  dt_masst = data.table(uid_leaf = c(), p_detected = c(), Cosine = c())
}



p_t = ggtree(tree, layout='circular', size=0.15, open.angle=5) + 
ggtitle(paste0(c('GNPS:', lib_id, 'Pubchem:', cid), collapse = ' ')) +
labs(caption = molecule_name) +
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
            mapping=aes(y=uid_leaf, color = Cosine)) +
            scale_color_viridis(option = "D", discrete = FALSE)


tree_data <- as.data.table(p_t$data)

# Join with tree_data to get x coordinates
dt_masst <- merge(dt_masst, tree_data, by.x = "uid_leaf", by.y = "label", all.x = TRUE)



  p_t = 
  p_t + geom_point2(data= dt_masst, aes(x = x, y = y, color = Cosine), alpha = 0.5) +
    scale_color_continuous(name = "Cosine")  +
    scale_color_viridis(option = "D", discrete = FALSE)
}

if ( nrow(dt_databases) > 0){

  dt_databases <- merge(dt_databases, tree_data, by.x = "uid_leaf", by.y = "label", all.x = TRUE)
  
  p_t = 
  p_t + geom_point2(data= dt_databases, aes(x = x, y = y, color = Database), color = 'red', shape = 4)

}

# Read the SVG content as a string
# svg_content <- readLines(args$mol_plot, warn = FALSE)

# if (!requireNamespace("ggsvg", quietly = TRUE)) {
#   remotes::install_github("coolbutuseless/ggsvg")
# }
# library(ggsvg)

# # Add the SVG to the plot
# p_t <- p_t  + annotation_custom(
#   grob = grid::grid.draw( svg_to_rasterGrob(svg_content)), 
#   xmin = Inf, xmax = Inf, 
#   ymin = Inf, ymax = Inf
# )


#svg_file <- image_read_svg(args$mol_plot)
#svg_file_scaled <- image_scale(svg_file, "3000x3000")
#svg_grob <- rasterGrob(as.raster(svg_file_scaled), interpolate = TRUE)
img <- image_read(args$mol_plot)
img_resized <- image_resize(img, "800x800")
img_raster <- rasterGrob(img, interpolate = TRUE)

# Overlay the image on the ggplot
p_t <- ggdraw(p_t) +
  draw_grob(img_raster, x = 0.95, y = 1, width = 0.4, height = 0.3, hjust = 1, vjust = 1)


# save_plot_with_svg(p_t, svg_grob, (paste0(c(lib_id, '_', cid, '.png'), collapse = '')))

ggsave(paste0(c('tree_',  lib_id, '_', cid, '.png'), collapse = ''), plot = p_t, width = 10, height = 8, dpi = 300)