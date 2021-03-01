# Phylogenetic tree plotting for https://github.com/pimbongaerts/pachyseris
#
# Author: Pim Bongaerts
# Copyright (C) 2021 Pim Bongaerts
# License: GPL

## Dependencies ==============================================================
suppressMessages(library("ggplot2"))
suppressMessages(library("ggtree"))
#suppressMessages(library("ggstance"))
#suppressMessages(library("reshape2"))
suppressMessages(library("dplyr"))
suppressMessages(library("grid"))
suppressMessages(library("gtable"))
suppressMessages(library("ape"))

## Main code =================================================================
### Functions ================================================================
plot_consensus_tree <- function(tree_data, popfile, offset){
  ggtree(tree_data) %<+% popfile +
    geom_tiplab(aes(color = factor(cluster)), align=TRUE, 
                size=3, linesize=.3, offset=offset) +
    hexpand(.1) +
    geom_tippoint(aes(color = factor(cluster)),  size=3) +
    scale_color_manual(values = color_scheme) +
    geom_point2(aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) > 95),
                color = "#000000") +
    geom_point2(aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) <= 95),
                shape = 4, size = 2) +
    theme(legend.position='none')
}

plot_density_tree <- function(tree_data, tip_order, popfile, offset){
  ggdensitree(tree_data, alpha = .3,
              layout = 'slanted', 
              tip.order = rev(tip_order)) %<+% popfile +
    geom_tippoint(aes(color = factor(cluster)),  size = 2) +
    scale_color_manual(values = color_scheme) +
    scale_x_reverse() +
    theme(legend.position='none')
}

plot_consensus_tree_compact <- function(tree_data, popfile, offset){
  ggtree(tree_data) %<+% popfile +
    geom_tiplab(aes(color = factor(cluster)), align=TRUE, 
                size=0, linesize=.3, offset=offset) +
    hexpand(.1) +
    scale_color_manual(values = color_scheme) +
    geom_point2(aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) > 95),
                shape = 16, color = "#000000") +
    geom_point2(aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label) <= 95),
                shape = 4, size = 2) +
    theme(legend.position='none')
}

tiporder <- function(tree_data){
  with(subset(fortify(tree_data), isTip), label[order(y, decreasing=T)])
}
### Initial set-up ===========================================================
# Read in datasets
setwd("/Users/pbongaerts/Dropbox/Data/pachyseris/server_data/B_RADseq/B6 - Phylogenomic analyses")

# Nexus tree
pachy_raxml_trees <- read.tree("raxml/RAxML_bootstrap.pachy_phylo_raxml")
pachy_raxml_cons_tree <- read.tree("raxml/RAxML_bipartitions.pachy_phylo_raxml")
pachy_raxml_cons_tree_rooted <- root(pachy_raxml_cons_tree, node = 51)

pachy_tetrad_trees <- read.tree("tetrad/pachy_tetrad_ref.tree.boots")
pachy_tetrad_cons_tree <- read.tree("tetrad/pachy_tetrad_ref.tree.cons")
pachy_tetrad_cons_tree_rooted <- root(pachy_tetrad_cons_tree, node = 79)
# Correct rerooting artefact:
# pachy_tetrad_cons_tree_rooted$node.label[2] = "100"

# Cluster assignments (from STRUCTURE K=6)
pachy_popfile = read.table("pachy_phylo_popfile2.txt",
                            col.names = c("sample", "cluster"), 
                            header = FALSE, stringsAsFactors = FALSE)

color_scheme = c("BL2" = "#08519c",  # BL2
                 "BLU" = "#377eb8",  # BLU
                 "GR2" = "#238b45",  # GR2 
                 "GRN" = "#4daf4a",  # GRN
                 "ISR" = "#810f7c",  # ISR
                 "RED" = "#e41a1c",  # RED
                 "RUG" = "black",
                 "INA" = "black",
                 "LEP" = "black",
                 "AGA" = "black")

### Produce individual plots =================================================

plot_raxml_cons_tree <-  plot_consensus_tree(pachy_raxml_cons_tree_rooted, pachy_popfile, 0)
plot_raxml_dens_tree <- plot_density_tree(pachy_raxml_trees,
                                          tiporder(pachy_raxml_cons_tree_rooted),
                                          pachy_popfile, 0)
plot_tetrad_cons_tree <- plot_consensus_tree(pachy_tetrad_cons_tree_rooted, pachy_popfile, 20)
plot_tetrad_dens_tree <- plot_density_tree(pachy_tetrad_trees, 
                                           tiporder(pachy_tetrad_cons_tree_rooted),
                                           pachy_popfile, 0)

plots_combined <- multiplot(plot_raxml_cons_tree, plot_tetrad_cons_tree,
                            plot_raxml_dens_tree, plot_tetrad_dens_tree, ncol = 2)

# Plots for supplementary
ggsave(plot_raxml_cons_tree, file="plot_raxml_cons_tree.pdf", width = 22, height = 22, units = "cm")
#ggsave(plot_raxml_dens_tree, file="plot_raxml_dens_tree.pdf", width = 22, height = 22, units = "cm")
ggsave(plot_tetrad_cons_tree, file="plot_tetrad_cons_tree.pdf", width = 22, height = 22, units = "cm")
ggsave(plot_tetrad_dens_tree, file="plot_tetrad_dens_tree.pdf", width = 22, height = 22, units = "cm")


# Plots for main figure
plot_raxml_cons_tree_compact <- plot_consensus_tree_compact(pachy_raxml_cons_tree_rooted, pachy_popfile, 0)
ggsave(plot_raxml_cons_tree_compact, file="plot_raxml_cons_tree_compact.pdf", width = 22, height = 11, units = "cm", useDingbats=FALSE)
ggsave(plot_raxml_dens_tree, file="plot_raxml_dens_tree_compact.pdf", width = 22, height = 11, units = "cm")
