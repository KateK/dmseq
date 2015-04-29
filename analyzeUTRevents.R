##----------------------------------
## ALE
##----------------------------------
gff_file <- "~/alt.splicing.annotations/hg19_Aug2014/ALE.hg19.gff3"

ALE_res <- analyzeUTRevents(gff_file)

write.table(ALE_res, file = "/home3/struck/Analysis_Projects/DMseq/data/ALE_metadata.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
## AFE
##----------------------------------
gff_file <- "~/alt.splicing.annotations/hg19_Aug2014/AFE.hg19.gff3"

AFE_res <- analyzeUTRevents(gff_file)

write.table(AFE_res, file = "/home3/struck/Analysis_Projects/DMseq/data/AFE_metadata.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
##----------------------------------
## Analysis Functions
##----------------------------------
##----------------------------------
analyzeUTRevents <- function(gff_file) {
    require(data.table)
    require(dplyr)
    require(GenomicRanges)
    
    gff <- fread(gff_file)

    gff <- gff %>% dplyr::select(V1, V4, V5, V9, V3, V7, V2)
    setNames(gff, c("chr", "start", "end", "isoform", "feature", "strand", "event_type"))

    gff_exons <- gff %>% dplyr::filter(feature == "exon")

    gff_exons[, isoform := gsub("ID=.+;", "", isoform)]
    gff_exons[, isoform := gsub("Parent=", "", isoform)]
    gff_exons[, gene_symbol := gsub("..$", "", isoform)]
    gff_exons[, exon_length := end - start + 1]

    gff_exons %<>% group_by(isoform) %>% mutate(isoform_length = sum(exon_length))

    gff_exons$seq <- appendSequences(gff_exons)
    
    return(gff_exons)
}

appendSequences <- function(bed) {
    require(BSgenome.Hsapiens.UCSC.hg19)
    sequences <- apply(bed, 1, function(x)
        {
            as.character(getSeq(Hsapiens, name = as.character(x[1]),
                                start = as.numeric(x[2]),
                                end = as.numeric(x[3]),
                                strand = as.character(x[6])))
        })    
    return(sequences)
}
