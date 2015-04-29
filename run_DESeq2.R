#!/usr/bin/env Rscript

##----------------------------------
## sample group
##----------------------------------
sample_group <- "tibialis"
# sample_group <- "quadricep"
# sample_group <- "heart"
# sample_group <- "bicep"

##----------------------------------
## base directories
##----------------------------------
base_dir <- paste("~/Analysis_Projects/DMseq/data/", sample_group, "/", sep = "")
counts_dir <- paste("~/alignment_counts/", sample_group, "/", sep = "")
output_dir <- paste("~/Analysis_Projects/DMseq/results/", sample_group, "/", sep = "")

##----------------------------------
## phenotypic data about samples
##----------------------------------
pdata_file <- paste(base_dir, paste("DM", sample_group, "pdata.txt", sep = "_"), sep = "")

##----------------------------------
## excluded samples
##----------------------------------
exclude_samples <- c("R153")

##----------------------------------
## read counts per gene 
##----------------------------------
countsTable_file <- paste(counts_dir, sample_group, "_countsTable.txt", sep = "")

##----------------------------------
## Load Data
##----------------------------------
require(DESeq2)
require(dplyr)

countsTable <- read.table(countsTable_file, header=TRUE, row.names=1, sep="\t")
countsTable <- countsTable[which(! rownames(countsTable) %in% c("alignment_not_unique", "ambiguous", "no_feature", "not_aligned", "too_low_aQual")),]
countsTable <- countsTable[, !names(countsTable) %in% exclude_samples]

pdata <- read.table(pdata_file, header=TRUE, row.names=1, sep="\t")
names(countsTable) <- rownames(pdata)
names(pdata)[1] <- "sample"
pdata <- pdata %>% filter(!sample %in% exclude_samples)

##----------------------------------
## Run DESeq2
## Grouping by Diagnosis
##----------------------------------
dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = pdata, design = ~ Diagnosis)
dds <- DESeq(dds)

## return normalized counts
# normalizedCountsTable <- as.data.frame(counts(dds, normalized=T))
# write.table(normalizedCountsTable, file="~/Documents/Research/DM_tibialis/normalizedCountsTable.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

plotMA(dds, pvalCutoff = 0.01)

res_DM1_vs_Normal <- as.data.frame(results(dds, contrast=c("Diagnosis", "DM1", "Normal")))
res_DM1_vs_Normal <- res_allDM1_vs_Normal[order(res_allDM1_vs_Normal$padj),]

output_file <- paste(output_dir, sample_group, "_DESeq2_DM1_vs_Normal.txt", sep = "")

write.table(res, file="~/Desktop/DESeqResults_DM1_vs_Normal.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

##----------------------------------
## Run DESeq2
## Grouping by strength
##----------------------------------
dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = pdata, design = ~ Strength_Grouping)
dds <- DESeq(dds)

res_W_vs_N <- as.data.frame(results(dds, contrast = c("Strength_Grouping", "Weak_DM1", "Control")))
res_W_vs_N <- res_W_vs_N[order(res_W_vs_N$padj),]

res_P_vs_N <- as.data.frame(results(dds, contrast = c("Strength_Grouping", "Proto_DM1", "Control")))
res_P_vs_N <- res_P_vs_N[order(res_P_vs_N$padj),]

res_NW_vs_N <- as.data.frame(results(dds, contrast = c("Strength_Grouping", "Not_Weak_DM1", "Control"))
res_NW_vs_N <- res_NW_vs_N[order(res_NW_vs_N$padj),]

# write.table(res_W_vs_N, file="~/Desktop/DESeqResults_W_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(res_NW_vs_N, file="~/Desktop/DESeqResults_NW_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(res_P_vs_N, file="~/Desktop/DESeqResults_P_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

