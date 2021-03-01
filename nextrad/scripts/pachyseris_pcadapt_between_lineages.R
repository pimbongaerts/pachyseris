library(pcadapt)
library(qvalue)
library(dplyr)

# Based on tutorial: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

# Comparisons
vcfs <- list("BLUvGRN" = read.pcadapt("pachy_pcadapt_BLU_vs_GRN.vcf", type = "vcf"),
             "BLUvRED" = read.pcadapt("pachy_pcadapt_BLU_vs_RED.vcf", type = "vcf"),
             "GRNvRED" = read.pcadapt("pachy_pcadapt_GRN_vs_RED.vcf", type = "vcf"),
             "BLUvISR" = read.pcadapt("pachy_pcadapt_BLU_vs_ISR.vcf", type = "vcf"),
             "GRNvISR" = read.pcadapt("pachy_pcadapt_GRN_vs_ISR.vcf", type = "vcf"),
             "REDvISR" = read.pcadapt("pachy_pcadapt_RED_vs_ISR.vcf", type = "vcf"))

popfiles <- list("BLUvGRN" = read.table("pachy_groups_BLU_GRN.txt")$V2,
                 "BLUvRED" = read.table("pachy_groups_BLU_RED.txt")$V2,
                 "GRNvRED" = read.table("pachy_groups_GRN_RED.txt")$V2,
                 "BLUvISR" = read.table("pachy_groups_BLU_ISR.txt")$V2,
                 "GRNvISR" = read.table("pachy_groups_GRN_ISR.txt")$V2,
                 "REDvISR" = read.table("pachy_groups_RED_ISR.txt")$V2)

# Conduct PCA using pcadapt and confirm that PC1  separates lineages only
# (and not e.g. regions within lineages)
pcadapt_outputs = list()
for (comparison in names(vcfs)){
  pcadapt_outputs[[as.name(comparison)]] <- pcadapt(input = vcfs[[as.name(comparison)]] , K = 2)
  plot(pcadapt_outputs[[as.name(comparison)]], option = "scores", pop = popfiles[[as.name(comparison)]])
}

# Rerun pcadapt with K = 1 and display results
for (comparison in names(vcfs)){
  pcadapt_outputs[[as.name(comparison)]] <- pcadapt(input = vcfs[[as.name(comparison)]] , K = 1)
  plot(pcadapt_outputs[[as.name(comparison)]], option = "scores", pop = popfiles[[as.name(comparison)]])
  plot(pcadapt_outputs[[as.name(comparison)]], option = "manhattan")
}

# Identify candidate SNPs
alpha <- 0.05 # expected false discovery rate lower than 5%
qvals = list()
for (comparison in names(vcfs)){
  qvals[[as.name(comparison)]] = qvalue(pcadapt_outputs[[as.name(comparison)]]$pvalues)$qvalues
  write.table(qvals[[as.name(comparison)]], file = paste(as.name(comparison), "_qvals.txt", sep = ""),
              row.names = TRUE, col.names = FALSE, quote = FALSE)
}