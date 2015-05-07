##----------------------------------
## ALE
##----------------------------------
gff_file <- "~/Projects/dmseq/data/event_annotations/ALE.hg19.gff3"

ALE_res <- UTRgffToBed(gff_file, event_type = "AFE")

write.table(ALE_res, file = "~/Projects/dmseq/data/ALE_metadata.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
## AFE
##----------------------------------
gff_file <- "~/Projects/dmseq/data/event_annotations/AFE.hg19.gff3"

AFE_res <- MisoGffToBed(gff_file, event_type = "ALE")

write.table(AFE_res, file = "~/Projects/dmseq/data/AFE_metadata.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
## Alt polyA
##----------------------------------
gff_file <- "~/Projects/dmseq/data/event_annotations/polyAseq.hg19.entrezID.v2.gff3"

polyA_res <- MisoGffToBed(gff_file, event_type = "polyA")

write.table(polyA_res, file = "~/Projects/dmseq/data/nonUTRevents.multi_metadata.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

##----------------------------------
## nonUTRevents
##----------------------------------
gff_file <- "~/Projects/dmseq/data/event_annotations/nonUTRevents.multi.hg19.gff3"

nonUTRevents_res <- MisoGffToBed(gff_file)

write.table(nonUTRevents_res, file = "~/Projects/dmseq/data/nonUTRevents.multi_metadata.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


##----------------------------------
## Analysis Functions
##----------------------------------
MisoGffToBed <- function(gff_file, event_type = "nonUTRevents") {
    require(data.table)
    require(dplyr)
    require(GenomicRanges)
    ## Read gff file and reformat to BED
    gff <- fread(gff_file)
    gff <- gff %>% dplyr::select(V1, V4, V5, V9, V3, V7, V2)
    setnames(gff, c("chr", "start", "end", "isoform", "feature", "strand", "event_type"))
    gff_exons <- gff %>% dplyr::filter(feature == "exon")
    ## Check event_type
    if (event_type == "nonUTRevents") {        
        gff_exons[, feature := gsub("ID=chr.+[+-].[A-Z].", "", gsub(";Parent=.+", "", isoform))]   
        gff_exons[, isoform := gsub("ID=.+;", "", isoform)]
        gff_exons[, isoform := gsub("Parent=", "", isoform)]
        gff_exons[, event_name := gsub("..$", "", isoform)]
        gff_exons[, exon_length := end - start + 1]
    } else {
        gff_exons[, isoform := gsub(";Parent=.+", "", isoform)]
        gff_exons[, isoform := gsub("ID=", "", isoform)]    
        gff_exons[, gene_symbol := gsub("..$", "", isoform)]
        gff_exons[, exon_length := end - start + 1]
        gff_exons <- gff_exons %>% group_by(isoform) %>% mutate(isoform_length = sum(exon_length))
    }     
    gff_exons$seq <- appendSequences(gff_exons)    
    return(gff_exons)
}

appendSequences <- function(bed) {
    require(BSgenome.Hsapiens.UCSC.hg19)
    sequences <- apply(bed, 1, function(x)
        {
            as.character(getSeq(Hsapiens,
                                name = as.character(x[1]),
                                start = as.numeric(x[2]),
                                end = as.numeric(x[3]),
                                strand = as.character(x[6])
                                )
                         )
        })    
    return(sequences)
}
