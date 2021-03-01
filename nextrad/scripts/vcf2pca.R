# PCA from VCF file
#
# Authors: Pim Bongaerts
# Copyright (C) 2019 Pim Bongaerts
# License: GPL
# Arguments: [vcf_filename] [number of axes to keep]

## Dependencies ==============================================================

library("adegenet")
library("vcfR")

## Main code =================================================================

## Obtain filenames through command-line arguments
args <- commandArgs(TRUE)
vcf_filename <- args[1]
axes <- args[2]

## Import with vcfR
coral.genind <- vcfR2genind(read.vcfR(vcf_filename))

## Set output files
filename_base <- strsplit(vcf_filename, "\\.")[[1]][[1]]
output_pca_filename <- paste(filename_base, "_pca.txt", sep = "")
output_pca_png_filename <- paste(filename_base, "_pca.png", sep = "")

## Principal component analysis
coral.scaled <- scaleGen(coral.genind, center=TRUE, scale=TRUE, NA.method="mean")  
pca1 <- dudi.pca(coral.scaled, cent = FALSE, scale = FALSE, scannf = FALSE,
                 nf = axes)
write.table(pca1$l1, file = output_pca_filename)

# Variance explained
eig.perc <- 100*pca1$eig/sum(pca1$eig)
eig.perc[1]
eig.perc[2]
eig.perc[3]