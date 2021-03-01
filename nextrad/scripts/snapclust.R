# Run snapclust for replicate VCF files and output to CSV
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

## Main code =================================================================
# Obtain filenames through command-line arguments
args <- commandArgs(TRUE)
arg_pattern <- args[1]   # e.g. "pachy_singlesnp*.*"
arg_cluster_start <- args[2]   # e.g. 4
arg_cluster_end <- args[3]     # e.g. 6
no_clusters = 1 + arg_cluster_end - arg_cluster_start

# Obtain list of vcf files to summarise
vcf_files <- list.files(path=".", pattern=arg_pattern, 
                        full.names=T, recursive=FALSE)

# Loop over files and run snapclust for each K
assignment_groups <- list()
assignment_probs <- list()
for(filename in vcf_files){
  vcf_genind <- vcfR2genind(read.vcfR(filename))
  filebase <- sub('\\..*$', '', basename(filename))
  replicate = substr(filebase, 16, 16) # specific to pachy_snp filename
  # Run snapclust for each k
  for(i in arg_cluster_start:arg_cluster_end){
    assignment = snapclust(vcf_genind, pop.ini = "ward", n.start = 50, 
                           n.start.kmeans = 50, k = i)
    run_name = paste(replicate, "_k", i, sep="")
    assignment_groups[[run_name]] = assignment$group
    assignment_probs[[run_name]] = assignment$proba
  }
}
# Combine individual runs into single data frame and save as csv
assignment_groups_combined <- as.data.frame(assignment_groups)
assignment_probs_combined <- as.data.frame(assignment_probs)
write.csv(assignment_groups_combined, file = 'snapclust_assignment_groups.csv')
write.csv(assignment_probs_combined, file = 'snapclust_assignment_probs.csv')

# Plots all assignments into a single figure
png(filename = "snapclust_summary.png", width = 2000, height = 1000)
assignment_groups_melt <- melt(as.matrix(assignment_groups_combined))
colnames(assignment_groups_melt) <- c("sample", "k", "assignment")
ggplot(assignment_groups_melt, aes(k, sample)) + 
  geom_tile(aes(fill = assignment)) +
  facet_grid(rows = vars(substr(sample, 1, 3)),
             scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 2))
dev.off()