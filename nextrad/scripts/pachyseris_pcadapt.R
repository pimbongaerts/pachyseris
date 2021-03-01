library(pcadapt)
library(qvalue)
library(dplyr)

# Based on tutorial: https://bcm-uga.github.io/pcadapt/articles/pcadapt.html

# Read VCF and popfiles
BLU_vcf <- read.pcadapt("pachy_BLU_b5a.vcf", type = "vcf")
GRN_vcf <- read.pcadapt("pachy_GRN_b5a.vcf", type = "vcf")
RED_vcf <- read.pcadapt("pachy_RED_b5a.vcf", type = "vcf")

BLU_pops <- substr(read.table("pachy_popfile_BLU_b5a.txt")$V2, 1, 2)
GRN_pops <- substr(read.table("pachy_popfile_GRN_b5a.txt")$V2, 1, 2)
RED_pops <- substr(read.table("pachy_popfile_RED_b5a.txt")$V2, 1, 2)

# Conduct PCA using pcadapt
pcadapt_BLU <- pcadapt(input = BLU_vcf, K = 20)
pcadapt_GRN <- pcadapt(input = GRN_vcf, K = 20)
pcadapt_RED <- pcadapt(input = RED_vcf, K = 20)

# Assess optimal choice for K through screeplot
plot(pcadapt_BLU, option = "screeplot") # K = 3 based on Cattell’s rule
plot(pcadapt_GRN, option = "screeplot") # K = 4 based on Cattell’s rule
plot(pcadapt_RED, option = "screeplot") # K = 2 based on Cattell’s rule

# Assess optimal choice for K through scoreplot
# BLU - no additional structure beyond K = 3
plot(pcadapt_BLU, option = "scores", pop = BLU_pops)
plot(pcadapt_BLU, option = "scores", i = 3, j = 4, pop = BLU_pops)
plot(pcadapt_BLU, option = "scores", i = 5, j = 6, pop = BLU_pops)
# GRN - no additional structure beyond K = 5 (so choosing K = 5)
plot(pcadapt_GRN, option = "scores", pop = GRN_pops)
plot(pcadapt_GRN, option = "scores", i = 3, j = 4, pop = GRN_pops)
plot(pcadapt_GRN, option = "scores", i = 5, j = 6, pop = GRN_pops)
# RED - no additional structure beyond K = 2
plot(pcadapt_RED, option = "scores", pop = RED_pops)
plot(pcadapt_RED, option = "scores", i = 3, j = 4, pop = RED_pops)
plot(pcadapt_RED, option = "scores", i = 5, j = 6, pop = RED_pops)

# Run PCA with final K
pcadapt_BLU <- pcadapt(input = BLU_vcf, K = 3)
pcadapt_GRN <- pcadapt(input = GRN_vcf, K = 5)
pcadapt_RED <- pcadapt(input = RED_vcf, K = 2)

# Display results
plot(pcadapt_BLU, option = "manhattan")
plot(pcadapt_GRN, option = "manhattan")
plot(pcadapt_RED, option = "manhattan")
plot(pcadapt_BLU, option = "qqplot")
plot(pcadapt_GRN, option = "qqplot")
plot(pcadapt_RED, option = "qqplot")
hist(pcadapt_BLU$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
hist(pcadapt_GRN$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
hist(pcadapt_RED$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")

# Candidate SNPs using an expected false discovery rate lower than 10%
alpha <- 0.1

BLU_qval <- qvalue(pcadapt_BLU$pvalues)$qvalues
write.table(BLU_qval, file="pachy_BLU_pcadapt_B5a.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE)
BLU_outliers <- which(BLU_qval < alpha)
count(get.pc(pcadapt_BLU, BLU_outliers), 'PC')

GRN_qval <- qvalue(pcadapt_GRN$pvalues)$qvalues
write.table(GRN_qval, file="pachy_GRN_pcadapt_B5a.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE)
GRN_outliers <- which(GRN_qval < alpha)
count(get.pc(pcadapt_GRN, GRN_outliers), 'PC')

RED_qval <- qvalue(pcadapt_RED$pvalues)$qvalues
write.table(RED_qval, file="pachy_RED_pcadapt_B5a.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE)
RED_outliers <- which(RED_qval < alpha)
count(get.pc(pcadapt_RED, RED_outliers), 'PC')

