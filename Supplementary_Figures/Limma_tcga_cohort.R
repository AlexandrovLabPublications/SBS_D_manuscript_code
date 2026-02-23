library(data.table)
library(ggplot2)
library(tidyverse)

# updated for different analysis levels
setwd('/Users/kazachkova/Documents/alexandrov_lab/mutographs/CRC_clonal_analysis/V3_Data/examine_SBS288D/Reviewers_Comments/Reviewers_comments_V2/TCGA_gene_expression_work/')


mrna_counts = read.csv("mrna_remove_dups_drop_lows.csv",row.names=1,header=T)
head(mrna_counts)


grouping = read.csv("metadata_tcga.csv",row.names=1,header=T)
grouping


in_both = intersect(colnames(mrna_counts),rownames(grouping))
length(in_both)
mrna_counts_ = subset(mrna_counts, select=in_both)

grouping <- grouping[in_both,]
grouping <- data.frame(grouping)
grouping$names <- in_both
rownames(grouping) <- grouping$names




# try limma 

library(limma)
library(edgeR) 
dge <- DGEList(counts = mrna_counts_)
grouping$SBS_D <- relevel(grouping$SBS_D, ref = "False")
design <- model.matrix(
  ~ gender + age_at_index + purity + subsite + SBS_D,
  data = grouping
)

colnames(design)
dge <- calcNormFactors(dge, method = "TMM")
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
colnames(fit$coefficients)
res_limma = topTable(
  fit,
  coef = "SBS_DTrue",   
  number = Inf,
  adjust.method = "BH"
)
write.csv(res_limma,'res_limma.csv')
