##----------------------------------
## required libraries
##----------------------------------
library(DESeq2)
library(dplyr)

##----------------------------------
## base directories
##----------------------------------
base_dir <- "~/Projects/DMseq/data/"
output_dir <- "~/Projects/DMseq/results/"

##----------------------------------
## phenotypic data about samples
##----------------------------------
pdata_file <- paste(base_dir, "DM_sample_pdata.txt", sep = "")

##----------------------------------
## excluded samples
##----------------------------------
exclude_samples <- c("R153")
# exclude_samples <- NULL

##----------------------------------
## read counts per gene 
##----------------------------------
countsTable_file <- paste(base_dir, "DMseq_allsamples_countsTable_geneSymbol.txt", sep = "")

##----------------------------------
## Load Data
##----------------------------------
countsTable <- read.table(countsTable_file, header=TRUE, row.names=1, sep="\t", check.names = FALSE)
countsTable <- countsTable[which(! rownames(countsTable) %in% c("__alignment_not_unique", "__ambiguous", "__no_feature", "__not_aligned", "__too_low_aQual")),]
names(countsTable) <- gsub("_counts_geneSymbol", "", names(countsTable))
countsTable <- countsTable[, !names(countsTable) %in% exclude_samples]

pdata <- read.table(pdata_file, header=TRUE, row.names=1, sep="\t")
pdata <- pdata[!rownames(pdata) %in% exclude_samples,]

##----------------------------------
## Subset data
##----------------------------------
sample_filter <- pdata$Tissue %in% c("Quad", "Tibialis") & rownames(pdata) %in% names(countsTable) & pdata$Diagnosis == "Control"
pdata <- pdata[sample_filter,]
countsTable <- countsTable[, names(countsTable) %in% rownames(pdata)]
countsTable <- countsTable[!is.na(rowSums(countsTable)),]
pdata <- pdata[names(countsTable),]

##----------------------------------
## Run DESeq2
## Grouping by Diagnosis
##----------------------------------
dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = pdata, design = ~Tissue) 
dds <- DESeq(dds)

## return normalized counts
# normalizedCountsTable <- as.data.frame(counts(dds, normalized=T))
# write.table(normalizedCountsTable, file="~/Documents/Research/DM_tibialis/normalizedCountsTable.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

plotMA(dds)

?plotMA

res <- as.data.frame(results(dds, contrast=c("Tissue", "Quad", "Tibialis")))
res <- res_quad_vs_tibialis[order(res_quad_vs_tibialis$padj),]

## res_DM1_vs_Normal <- as.data.frame(results(dds, contrast=c("Diagnosis", "DM1", "Normal")))
## res_DM1_vs_Normal <- res_allDM1_vs_Normal[order(res_allDM1_vs_Normal$padj),]

output_file <- paste(output_dir, "DESeq2_Control_Quad_vs_Tibialis.txt", sep = "")

write.table(res, file=output_file, sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

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

res_NW_vs_N <- as.data.frame(results(dds, contrast = c("Strength_Grouping", "Not_Weak_DM1", "Control")))
res_NW_vs_N <- res_NW_vs_N[order(res_NW_vs_N$padj),]

# write.table(res_W_vs_N, file="~/Desktop/DESeqResults_W_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(res_NW_vs_N, file="~/Desktop/DESeqResults_NW_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(res_P_vs_N, file="~/Desktop/DESeqResults_P_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

##----------------------------------
## Gene Ontology for significant results
##----------------------------------
library(GO.db)
library(org.Hs.eg.db)

# sig.res.pdata <- select(org.Hs.eg.db, keys = rownames(DE.AS[which(DE.AS$padj <= 0.05),]), columns=c("GENENAME", "GO", "ONTOLOGY"), keytype="SYMBOL")
sig.res.pdata <- select(org.Hs.eg.db, keys = rownames(res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1),]), columns=c("GENENAME", "GO", "ONTOLOGY"), keytype="SYMBOL")
sig.res.pdata$GO[is.na(sig.res.pdata$GO)] <- "GO:NA"
sig.res.pdata$ONTOLOGYTERM <- select(GO.db, keys = sig.res.pdata$GO, columns = c("TERM"), keytype = "GOID")$TERM

textplot(subset(sig.res.pdata, grepl("splicing", ONTOLOGYTERM)))
