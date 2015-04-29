#!/usr/bin/env Rscript

##----------------------------------
## load analysis functions
##----------------------------------
source("~/Analysis_Projects/DMseq/bin/summarizeMISOoutput_dev.R")

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
## Run DESeq2
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
                           
## By Diagnosis only
dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = pdata, design = ~ Diagnosis)

dds <- DESeq(dds)
# normalizedCountsTable <- as.data.frame(counts(dds, normalized=T))
# write.table(normalizedCountsTable, file="~/Documents/Research/DM_tibialis/normalizedCountsTable.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

plotMA(dds, pvalCutoff = 0.01)

res_allDM1_vs_Normal <- as.data.frame(results(dds, contrast=c("Diagnosis", "DM1", "Normal")))
res_allDM1_vs_Normal <- res_allDM1_vs_Normal[order(res_allDM1_vs_Normal$padj),]

# write.table(res, file="~/Desktop/DESeqResults_allDM1_vs_Normal.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

### By Diagnosis and Strength
pdata$Strength_Grouping  <- as.factor(unlist(lapply(pdata$Actual_Strength_6pt_Scale_Norm7_proto6.5, function(x) {
  if(x < 5) { return("W") } else { if(x %in% c(5,6)) { return("NW") } else { if(x == 6.5) { return("P") } else{ if(x == 7) { return("N")} }}}
})))

dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = pdata, design = ~ Strength_Grouping)
dds <- DESeq(dds, betaPrior=FALSE)

res_W_vs_N <- as.data.frame(results(dds, contrast = c("Strength_Grouping", "W", "N")))
res_W_vs_N <- res_W_vs_N[order(res_W_vs_N$padj),]

res_P_vs_N <- as.data.frame(results(dds, contrast = c("Strength_Grouping", "P", "N")))
res_P_vs_N <- res_P_vs_N[order(res_P_vs_N$padj),]

res_NW_vs_N <- as.data.frame(results(dds, contrast = c("Strength_Grouping", "NW", "N"))
res_NW_vs_N <- res_NW_vs_N[order(res_NW_vs_N$padj),]

# write.table(res_W_vs_N, file="~/Desktop/DESeqResults_W_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(res_NW_vs_N, file="~/Desktop/DESeqResults_NW_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)
# write.table(res_P_vs_N, file="~/Desktop/DESeqResults_P_vs_N.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

### Calc RPKMS
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(edgeR)
library(data.table)

# Load the annotation and reduce it
GTFfile <- "~/Projects/annotations_hg19/Homo_sapiens.GRCh37.75.gtf"
GTF <- import(GTFfile, format="gtf", asRangedData=F, feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_name))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_name <- rep(names(grl), elementLengths(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

DT <- data.table(as.data.frame(reducedGTF, row.names=NULL), key="gene_name")
reducedDT <- DT[, list(width=sum(width)), by = gene_name]
reducedDT <- reducedDT[which(reducedDT$gene_name %in% rownames(countsTable)),]

rpkms <- rpkm(countsTable[which(rownames(countsTable) %in% reducedDT$gene_name),], gene.length = reducedDT$width, normalized.lib.sizes=FALSE)

write.table(rpkms, file="~/Documents/Research/DM_tibialis/rpkmTable.txt", sep= "\t", row.names = TRUE, col.names = TRUE, quote=FALSE)

### Other Analyses
sig.res <- res[which(res$padj < 0.05),]
up <- res[which(res$log2FoldChange >= 1 & res$padj < 0.05),]
down <- res[which(res$log2FoldChange <= -1 & res$padj < 0.05),]

plotMA(dds, pvalCutoff = 0.05, ylim=c(-4,4))

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
plotPCA(vsd, intgroup="Diagnosis")

rld <- rlogTransformation(dds, blind=TRUE)
distsRL <- dist(t(assay(rld)))

library("RColorBrewer")
library("gplots")

hc <- hclust(distsRL, method = "ward", members = NULL)
plot(hc)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- paste(rownames(colData(dds)), colData(dds)$Diagnosis, sep=" : ")
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(6, 6))

countsTable.sums <- rowSums(countsTable)
notAllZero <- (rowSums(counts(dds))>0)
counts.subset <- countsTable[which(countsTable.sums >= 500),]
count.cors <- cor(counts.subset, method="pearson", use="na.or.complete")
pheatmap(count.cors)

### Gene Ontology for significant results

library(GO.db)
library(org.Hs.eg.db)

# sig.res.pdata <- select(org.Hs.eg.db, keys = rownames(DE.AS[which(DE.AS$padj <= 0.05),]), columns=c("GENENAME", "GO", "ONTOLOGY"), keytype="SYMBOL")
sig.res.pdata <- select(org.Hs.eg.db, keys = rownames(res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1),]), columns=c("GENENAME", "GO", "ONTOLOGY"), keytype="SYMBOL")
sig.res.pdata$GO[is.na(sig.res.pdata$GO)] <- "GO:NA"
sig.res.pdata$ONTOLOGYTERM <- select(GO.db, keys = sig.res.pdata$GO, columns = c("TERM"), keytype = "GOID")$TERM
subset(sig.res.pdata, grepl("RNA splicing", ONTOLOGYTERM))
