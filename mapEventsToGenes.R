##-------------------------------
## required libraries
##-------------------------------
library(data.table)
library(GenomicRanges)

##-------------------------------
## map splicing events to genes
##-------------------------------
gtf <- fread("~/Projects/annotations_hg19/Homo_sapiens.GRCh37.75.gtf")

genes <- gtf[gtf$V3 == "gene",]
genes$V6 <- unlist(lapply(strsplit(genes$V9, "[; \"]"), function(x) x[8]))
genes <- genes %>% dplyr::select(V1, V4, V5, V6, V8, V7)
setNames(genes, c("chr", "start", "end", "name", "score", "strand"))

chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                 "chr22", "chrM", "chrX", "chrY")

genes <- genes[chr %in% chromosomes]

nonUTRevents <- fread("~/Projects/DMseqAnalysis/data/event_annotations/nonUTRevents.multi.hg19.gff3")
nonUTRevents <- nonUTRevents %>% dplyr::select(V1, V4, V5, V9, V3, V7, V2)
setNames(nonUTRevents, c("chr", "start", "end", "name", "score", "strand", "class"))

nonUTRevents_exons <- nonUTRevents %>% dplyr::filter(score == "exon")
nonUTRevents_exons$name <- gsub("ID=.+;", '', nonUTRevents_exons$name)
nonUTRevents_exons$name <- gsub("Parent=", "", nonUTRevents_exons$name)
nonUTRevents_exons$name <- gsub("Name=", "", nonUTRevents_exons$name)
nonUTRevents_exons$name <- gsub("..$", "", nonUTRevents_exons$name)
nonUTRevents_exons <- unique(nonUTRevents_exons)

genes_granges <- makeGRangesFromDataFrame(as.data.frame(genes), keep.extra.columns = TRUE)
nonUTRevents_granges <- makeGRangesFromDataFrame(as.data.frame(nonUTRevents_exons), keep.extra.columns = TRUE)

hits <- findOverlaps(nonUTRevents_granges, genes_granges)

eventKey <- data.frame(event_name = as.data.frame(nonUTRevents_granges[queryHits(hits),])$name,
                       gene_symbol = as.data.frame(genes_granges[subjectHits(hits),])$name)
eventKey <- unique(eventKey)

eventKey <- eventKey %>% group_by(event_name) %>% summarize(name = paste(gene_symbol, collapse = ","))

write.table(eventKey, file = "nonUTRevents.multi.hg19.to.geneSymbol.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
