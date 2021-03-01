# Evaluate clusters with snapclust for replicate VCF files
#
#
# Authors: Pim Bongaerts
# Copyright (C) 2019 Pim Bongaerts
# License: GPL
# Arguments: [vcf_filepattern]

## Dependencies ==============================================================

suppressMessages(library("adegenet"))
suppressMessages(library("vcfR"))
suppressMessages(library("ggplot2"))
suppressMessages(library("reshape2"))
suppressMessages(library("gridExtra"))

## Functions =================================================================
ic_plot <- function(ic_data, ic_label){
  ggplot(ic_data, aes(K,IC)) +
    geom_point(aes(colour = rep)) +
    geom_line(aes(colour = rep)) +
    theme_bw() +
    ylab(ic_label) +
    theme(legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size=12),
          axis.title = element_text(size=12),
          axis.line = element_line(colour = "black"))
}

## Main code =================================================================

# Obtain filenames through command-line arguments
args <- commandArgs(TRUE)
arg_pattern <- args[1]  #e.g. "pachy_singlesnp*.*""

## Set output file
output_png_filename <- "snapclust_pachy.png"

# Obtain list of vcf files to summarise
vcf_files <- list.files(path=".", pattern=arg_pattern, 
                        full.names=T, recursive=FALSE)

# Loop over files and add output data to list
aic_results <- list()
bic_results <- list()
kic_results <- list()
for(filename in vcf_files){
  vcf_genind<- vcfR2genind(read.vcfR(filename))
  aic_results[[filename]] <- snapclust.choose.k(20, vcf_genind, IC = AIC)
  bic_results[[filename]] <- snapclust.choose.k(20, vcf_genind, IC = BIC)
  kic_results[[filename]] <- snapclust.choose.k(20, vcf_genind, IC = KIC)
}

# Plot results
aic_melt <- melt(as.matrix(as.data.frame(aic_results)))
colnames(aic_melt) <- c("K", "rep", "IC")
write.table(aic_melt, file = "pachy_aic.txt")
aic_plot <- ic_plot(aic_melt, "AIC")

bic_melt <- melt(as.matrix(as.data.frame(bic_results)))
colnames(bic_melt) <- c("K", "rep", "IC")
write.table(bic_melt, file = "pachy_bic.txt")
bic_plot <- ic_plot(bic_melt, "BIC")

kic_melt <- melt(as.matrix(as.data.frame(kic_results)))
colnames(kic_melt) <- c("K", "rep", "IC")
write.table(kic_melt, file = "pachy_kic.txt")
kic_plot <- ic_plot(kic_melt, "KIC")

png(filename = output_png_filename, width = 1300, height = 350)
grid.arrange(aic_plot, bic_plot, kic_plot, nrow = 1, ncol=3)
dev.off()
