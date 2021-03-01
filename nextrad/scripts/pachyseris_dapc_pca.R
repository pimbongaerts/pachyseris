# DAPC and PCA for Pachyseris dataset
#
# Authors: Pim Bongaerts
# Copyright (C) 2020 Pim Bongaerts
# License: GPL

## Dependencies ==============================================================
suppressMessages(library("adegenet"))
suppressMessages(library("vcfR"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("reshape2"))
#suppressMessages(library("gridExtra"))
suppressMessages(library("cowplot"))

## Functions =================================================================
conduct_pca <- function(genind_data){
  genind_scaled <- scaleGen(genind_data, center=TRUE, 
                            scale=TRUE, NA.method="mean")  
  dudi.pca(genind_scaled, cent = FALSE, scale = FALSE, 
           scannf = FALSE, nf = 2)
}

conduct_dapc <- function(genind_data){
  max_PCAs <- as.integer(length(genind_data$pop) / 3)
  dapc_pop <- dapc(genind_data, n.pca = max_PCAs, n.da = 6)
  optimum_score <- optim.a.score(dapc_pop)
  dapc(genind_data, n.pca = optimum_score$best, n.da = 6)
}

set_pop_to_subregion <- function(pachy_genind){
  new_genind <- pachy_genind
  pop(new_genind) <- substr(rownames(new_genind@tab), 3, 4)
  pop(new_genind) <- recode_pop2region(pop(new_genind))
  return(new_genind)
}

set_pop_to_depth <- function(pachy_genind){
  new_genind <- pachy_genind
  pop(new_genind) <- paste(substr(rownames(new_genind@tab), 3, 3), 
                           substr(rownames(new_genind@tab), 5, 5), sep = "")
  return(new_genind)
}

remove_samples <- function(pachy_genind, groups_to_remove){
  pachy_genind[i = !grepl(groups_to_remove, rownames(pachy_genind@tab))]
}

recode_pop2region <- function(pops){
  recode(pops, "PG" = "PN", "PV"= "PN", 
         "GS" = "GF", "GG" = "GF", "GT" = "GF", "GY" = "GF", 
         "GD" = "GN", "GR" = "GN", "GU" = "GF", "GA" = "GF", "GM" = "GM", 
         "CF" = "CF", "CH" = "CH", "CB" = "CB", "CO" = "CO", "CY" = "CO")
}

## Plot functions ============================================================
compo_plot <- function(posterior_data, colors){
  # Organize data for plotting
  compo_data <- as.data.frame.matrix(posterior_data)
  compo_data$sample <- rownames(compo_data)
  compo_data$pop <- recode_pop2region(substr(compo_data$sample, 3, 4))
  compo_data$pop <- factor(compo_data$pop, levels = pop_order)
  compo_data <- melt(compo_data, id.vars = c("sample", "pop"))
  # Plot stacked bar graphs
  ggplot(compo_data, aes(x = sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "fill", width = 1, 
             size = 0.25, colour = "black") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(. ~ pop, space = "free", 
               scales = "free", switch = "x") +
    guides(fill = FALSE) +
    xlab("") +
    ylab("") +
    theme(plot.title = element_text(size = 8, 
                                    vjust = 0, 
                                    face = "bold"),
          #axis.text.x = element_blank(), 
          axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          panel.margin.x = unit(0.25, "lines"),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      size = 1, linetype="solid"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 8, angle = 270))
}

depth_compo_plot <- function(posterior_data, colors){
  # Organize data for plotting
  compo_data <- as.data.frame.matrix(posterior_data)
  compo_data$sample <- rownames(compo_data)
  compo_data$pop <- paste(substr(compo_data$sample, 3, 3),
                                 substr(compo_data$sample, 5, 5), sep = "")
  compo_data <- melt(compo_data, id.vars = c("sample", "pop"))
  # Plot stacked bar graphs
  ggplot(compo_data, aes(x = sample, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "fill", width = 1, 
             size = 0.25, colour = "black") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(. ~ pop, space = "free", 
               scales = "free", switch = "x") +
    guides(fill = FALSE) +
    xlab("") +
    ylab("") +
    theme(plot.title = element_text(size = 8, 
                                    vjust = 0, 
                                    face = "bold"),
          axis.text.x = element_blank(), 
          #axis.text.x = element_text(size = 6, angle = 90, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_blank(),
          panel.margin.x = unit(0.25, "lines"),
          panel.border = element_rect(fill = NA, colour = "black", 
                                      size = 1, linetype="solid"),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.text.y = element_text(size = 8, angle = 270))
}

pca_plot <- function(pca_data, colors){
  # Organize data for plotting
  pca_data$region <- substr(rownames(pca_data), 3, 3)
  pca_data$region <- recode(pca_data$region, "H" = "G")
  centroids <- aggregate(cbind(RS1, RS2) ~ region, data = pca_data, mean)
  pca_data <- merge(pca_data, centroids, by="region", suffixes=c("",".centroid"))
  # Scatter plot with connected centroid clusters
  ggplot(pca_data) +
    geom_point(aes(x=RS1, y=RS2, color=region), size = 1) +
    geom_point(data=centroids, aes(x=RS1, y=RS2, color=region), size = 2) +
    geom_segment(aes(x=RS1.centroid, y=RS2.centroid, xend=RS1, yend=RS2, color=region), size = 0.2) +
    scale_colour_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(plot.title = element_text(size = 8, vjust = 0, face="bold"),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

dapc_plot <- function(indcoor_data, colors){
  # Organize data for plotting
  dapc_data <- as.data.frame.matrix(indcoor_data)
  dapc_data$sample <- rownames(dapc_data)
  dapc_data$regiondepth <- paste(substr(dapc_data$sample, 3, 3),
                                 substr(dapc_data$sample, 5, 5), sep = "")
  
  centroids <- aggregate(cbind(LD1, LD2) ~ regiondepth, data = dapc_data, mean)
  dapc_data <- merge(dapc_data, centroids, by="regiondepth", suffixes=c("",".centroid"))
  # Plot stacked bar graphs  
  ggplot(dapc_data) +
    geom_point(aes(x=LD1, y=LD2, color=regiondepth), size = 1) +
    geom_point(data=centroids, aes(x=LD1, y=LD2, color=regiondepth), size = 2) +
    geom_segment(aes(x=LD1.centroid, y=LD2.centroid, xend=LD1, yend=LD2, color=regiondepth), size = 0.2) +
    scale_colour_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_bw() +
    theme(plot.title = element_text(size = 8, vjust = 0, face="bold"),
          legend.title = element_blank(),
          legend.position = "none",
          axis.text = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

dapc_loadingplot <- function(dapc_loadings, loading_color){
  # Organize data for plotting
  dapc_loadings_melt <- melt(dapc_loadings)
  colnames(dapc_loadings_melt) <- c("locus", "axis", "loading")
  dapc_loadings_melt$loading <- with(dapc_loadings_melt, 
                                     ifelse(axis == "LD1", loading, -loading))
  
  # Bar graph with DAPC loadings
  ggplot(dapc_loadings_melt, aes(x = locus, y = loading)) +  
    geom_bar(stat="identity", fill = loading_color) +
    facet_grid(axis ~ ., space = "free", scales = "free") +
    xlab("") +
    ylab("Loadings") +
    #scale_fill_manual(values = colours) +
    theme(legend.position="none",
          axis.text.x =  element_blank(),
          axis.title.x =  element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          #panel.margin.x = unit(0.2, "lines"),
          panel.margin.y = unit(0, "lines"),
    )
}
## Main code =================================================================

### Plot settings  ===========================================================
reanalyze = FALSE
pop_order <- rev(c("PN", "GF", "GN", "GM", "CF", "CH", "CB", "CO"))

red_pca_colors <- c("C" = "#ae017e", "G" = "#cb181d")
grn_pca_colors <- c("C" = "#02818a", "G" = "#238b45", "P" = "#a1dab4")
blu_pca_colors <- c("C" = "#6a51a3", "G" = "#2171b5", "P" = "#41b6c4")

red_compo_colors <- c("CO" = "#ae017e", "CB" = "#f768a1", "CF" = "#fbb4b9",
                      "GF" = "#cb181d", "GN" = "#fb6a4a", "GM" = "#fcae91")
grn_compo_colors <- c("CO" = "#02818a", "CB" = "#67a9cf", "CF" = "#bdc9e1",
                      "GF" = "#238b45", "GN" = "#74c476", "GM" = "#bae4b3",
                      "PN" = "#a1dab4")
blu_compo_colors <- c("CO" = "#6a51a3", "CB" = "#9e9ac8",
                      "GF" = "#2171b5", "GN" = "#6baed6", "GM" = "#3182bd",
                      "PN" = "#41b6c4")

red_dapc_colors <- c("GB" = "#fcae91", "GS" = "#fb6a4a", "GD" = "#cb181d")
grn_dapc_colors <- c("GB" = "#bae4b3", "GS" = "#74c476", "GD" = "#238b45")
blu_dapc_colors <- c("GB" = "#bdd7e7", "GS" = "#6baed6", "GD" = "#2171b5")

# Define matrix with overall layout (x, y, width, height)
screen_coords <- rbind(c( 0/18,  16/24, 10/18, 8/24),  # grn_compo
                       c( 0/18,  8/24, 10/18, 8/24),  # red_compo
                       c( 0/18,  0/24, 10/18, 8/24),  # blu_compo
                       c( 10/18, 16/24, 4/18,  8/24),  # grn_pca
                       c( 10/18, 8/24, 4/18,  8/24),  # red_pca
                       c( 10/18, 0/24, 4/18,  8/24),  # blu_pca
                       c( 14/18, 16/24, 4/18,  8/24),  # grn_density
                       c( 14/18, 8/24, 4/18,  8/24),  # red_density
                       c( 14/18, 0/24, 4/18,  8/24))  # blu_density

colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("grn_compo", "red_compo", "blu_compo", 
                             "grn_pca", "red_pca", "blu_pca",
                             "grn_dapc", "red_dapc", "blu_dapc")

### Read datasets ============================================================
vcf_red = "pachy_RED_b5a.vcf"
vcf_grn = "pachy_GRN_b5a.vcf"
vcf_blu = "pachy_BLU_b5a.vcf"

if (reanalyze == TRUE) {
  # red datasets
  red_genind <- vcfR2genind(read.vcfR(vcf_red))
  red_genind_subregion <- set_pop_to_subregion(red_genind)
  # TODO: removing Holmes Reef (C) and PNG as very few samples & 2 PCA outliers
  red_genind_subregion <- remove_samples(red_genind_subregion, "PSCH|PSP|PSCYSH8356|PSGABX7659")
  red_genind_depth <- set_pop_to_depth(red_genind)
  red_genind_depth <- remove_samples(red_genind_depth, "PSC|PSP")
  #table(pop(red_genind_depth))
  
  # grn datasets
  grn_genind <- vcfR2genind(read.vcfR(vcf_grn))
  grn_genind_subregion <- set_pop_to_subregion(grn_genind)
  # TODOP: removing Holmes Reef (C) and sperm sample:
  grn_genind_subregion <- remove_samples(grn_genind_subregion, "PSCH|PSHS")
  grn_genind_depth <- set_pop_to_depth(grn_genind)
  grn_genind_depth <- remove_samples(grn_genind_depth, "PSHS|PSC|PSP")
  #table(pop(grn_genind_depth))
  
  # blu datasets
  blu_genind <- vcfR2genind(read.vcfR(vcf_blu))
  blu_genind_subregion <- set_pop_to_subregion(blu_genind)
  # TODOP: removing Holmes Reef (C) and sperm sample:
  blu_genind_subregion <- remove_samples(blu_genind_subregion, "PSCF|PSHS")
  blu_genind_depth <- set_pop_to_depth(blu_genind)
  blu_genind_depth <- remove_samples(blu_genind_depth, "PSP|PSC|PSP")
  #table(pop(blu_genind_depth))
  
  ### Conduct structuring analyses =============================================
  red_pca <- conduct_pca(red_genind_subregion)
  red_dapc_subregion <- conduct_dapc(red_genind_subregion)
  red_dapc_depth<- conduct_dapc(red_genind_depth)
  
  grn_pca <- conduct_pca(grn_genind_subregion)
  grn_dapc_subregion <- conduct_dapc(grn_genind_subregion)
  grn_dapc_depth <- conduct_dapc(grn_genind_depth)
  
  blu_pca <- conduct_pca(blu_genind_subregion)
  blu_dapc_subregion <- conduct_dapc(blu_genind_subregion)
  blu_dapc_depth <- conduct_dapc(blu_genind_depth)
}
### Produce plots ============================================================
red_compo_plot <- compo_plot(red_dapc_subregion$posterior, red_compo_colors)
red_pca_plot <- pca_plot(red_pca$l1, red_pca_colors)
red_dapc_plot <- dapc_plot(red_dapc_depth$ind.coord, red_dapc_colors)

grn_compo_plot <- compo_plot(grn_dapc_subregion$posterior, grn_compo_colors)
grn_pca_plot <- pca_plot(grn_pca$l1, grn_pca_colors)
grn_dapc_plot <- dapc_plot(grn_dapc_depth$ind.coord, grn_dapc_colors)

blu_compo_plot <- compo_plot(blu_dapc_subregion$posterior, blu_compo_colors)
blu_pca_plot <- pca_plot(blu_pca$l1, blu_pca_colors)
blu_dapc_plot <- dapc_plot(blu_dapc_depth$ind.coord, blu_dapc_colors)

### Composite of individual plots ============================================

ggdraw() +
  draw_plot(blu_compo_plot, 
            screen_coords["blu_compo", "x"], screen_coords["blu_compo", "y"], 
            screen_coords["blu_compo", "width"], screen_coords["blu_compo", "height"]) +
  draw_plot(red_compo_plot, 
            screen_coords["red_compo", "x"], screen_coords["red_compo", "y"], 
            screen_coords["red_compo", "width"], screen_coords["red_compo", "height"]) +
  draw_plot(grn_compo_plot, 
            screen_coords["grn_compo", "x"], screen_coords["grn_compo", "y"], 
            screen_coords["grn_compo", "width"], screen_coords["grn_compo", "height"]) +
  draw_plot(grn_pca_plot, 
            screen_coords["grn_pca", "x"], screen_coords["grn_pca", "y"], 
            screen_coords["grn_pca", "width"], screen_coords["blu_pca", "height"]) +
  draw_plot(red_pca_plot, 
            screen_coords["red_pca", "x"], screen_coords["red_pca", "y"], 
            screen_coords["red_pca", "width"], screen_coords["red_pca", "height"]) +
  draw_plot(blu_pca_plot, 
            screen_coords["blu_pca", "x"], screen_coords["blu_pca", "y"], 
            screen_coords["blu_pca", "width"],screen_coords["blu_pca", "height"]) +
  draw_plot(grn_dapc_plot, 
            screen_coords["grn_dapc", "x"], screen_coords["grn_dapc", "y"], 
            screen_coords["grn_dapc", "width"], screen_coords["grn_dapc", "height"]) +
  draw_plot(red_dapc_plot, 
            screen_coords["red_dapc", "x"], screen_coords["red_dapc", "y"], 
            screen_coords["red_dapc", "width"], screen_coords["red_dapc", "height"]) +
  draw_plot(blu_dapc_plot, 
            screen_coords["blu_dapc", "x"], screen_coords["blu_dapc", "y"], 
            screen_coords["blu_dapc", "width"], screen_coords["blu_dapc", "height"])

ggsave("figure.pdf", width = 24, height = 16.8, units = "cm")
ggsave("figure.png", width = 24, height = 16.8, units = "cm")

### Additional supplementary plot ============================================
screen_coords <- rbind(c( 0/18,  16/24, 6/18, 8/24),  # grn_compo
                       c( 0/18,  8/24,  6/18, 8/24),  # red_compo
                       c( 0/18,  0/24,  6/18, 8/24),  # blu_compo
                       c( 6/18, 16/24,  4/18,  8/24),  # grn_dapc
                       c( 6/18, 8/24,   4/18,  8/24),  # red_dapc
                       c( 6/18, 0/24,   4/18,  8/24),  # blu_dapc
                       c( 10/18, 16/24, 8/18,  8/24),  # grn_l1
                       c( 10/18, 8/24,  8/18,  8/24),  # red_l1
                       c( 10/18, 0/24,  8/18,  8/24),  # blu_l1
                       c( 14/18, 16/24, 4/18,  8/24),  # grn_l2
                       c( 14/18, 8/24,  4/18,  8/24),  # red_l2
                       c( 14/18, 0/24,  4/18,  8/24))  # blu_l2
colnames(screen_coords) <- c("x", "y", "width", "height")
rownames(screen_coords) <- c("grn_compo", "red_compo", "blu_compo",
                             "grn_dapc", "red_dapc", "blu_dapc",
                             "grn_l1", "red_l1", "blu_l1",
                             "grn_l2", "red_l2", "blu_l2")
ggdraw() +
  draw_plot(depth_compo_plot(blu_dapc_depth$posterior, blu_dapc_colors), 
            screen_coords["blu_compo", "x"], screen_coords["blu_compo", "y"], 
            screen_coords["blu_compo", "width"], screen_coords["blu_compo", "height"]) +
  draw_plot(depth_compo_plot(red_dapc_depth$posterior, red_dapc_colors), 
            screen_coords["red_compo", "x"], screen_coords["red_compo", "y"], 
            screen_coords["red_compo", "width"], screen_coords["red_compo", "height"]) +
  draw_plot(depth_compo_plot(grn_dapc_depth$posterior, grn_dapc_colors), 
            screen_coords["grn_compo", "x"], screen_coords["grn_compo", "y"], 
            screen_coords["grn_compo", "width"], screen_coords["grn_compo", "height"]) +
  draw_plot(grn_dapc_plot, 
            screen_coords["grn_dapc", "x"], screen_coords["grn_dapc", "y"], 
            screen_coords["grn_dapc", "width"], screen_coords["grn_dapc", "height"]) +
  draw_plot(red_dapc_plot, 
            screen_coords["red_dapc", "x"], screen_coords["red_dapc", "y"], 
            screen_coords["red_dapc", "width"], screen_coords["red_dapc", "height"]) +
  draw_plot(blu_dapc_plot, 
            screen_coords["blu_dapc", "x"], screen_coords["blu_dapc", "y"], 
            screen_coords["blu_dapc", "width"],screen_coords["blu_dapc", "height"]) +
  draw_plot(dapc_loadingplot(grn_dapc_depth$var.contr, "#238b45"),
            screen_coords["grn_l1", "x"], screen_coords["grn_l1", "y"], 
            screen_coords["grn_l1", "width"], screen_coords["grn_l1", "height"]) +
  draw_plot(dapc_loadingplot(red_dapc_depth$var.contr, "#cb181d"),
            screen_coords["red_l1", "x"], screen_coords["red_l1", "y"], 
            screen_coords["red_l1", "width"], screen_coords["red_l1", "height"]) +
  draw_plot(dapc_loadingplot(blu_dapc_depth$var.contr, "#2171b5"),
            screen_coords["blu_l1", "x"], screen_coords["blu_l1", "y"], 
            screen_coords["blu_l1", "width"], screen_coords["blu_l1", "height"])

ggsave("supp_figure.pdf", width = 24, height = 16.8, units = "cm")
ggsave("supp_figure.png", width = 24, height = 16.8, units = "cm")
