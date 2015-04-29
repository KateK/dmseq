##----------------------------------
## required libraries
##----------------------------------
library(data.table)
library(dplyr)
library(tidyr)
library(GenomicRanges)

##----------------------------------
## load data
##----------------------------------
miso_bed <- tbl_dt(fread("~/Projects/DMseq/data/event_annotations/nonUTRevents.multi.hg19.bed"))
array_events <- tbl_dt(fread("~/Projects/DMseq/data/tibialis/arrayevents_bed.txt"))

## adjust events where the start > end 
tmp <- miso_bed %>% filter(start > end) %>% select(chr, end, start, name, score, strand)
setNames(tmp, c("chr", "start", "end", "name", "score", "strand"))
miso_bed <- tbl_dt(rbind(tmp, miso_bed %>% filter( end > start)))

array_events_granges <- makeGRangesFromDataFrame(as.data.frame(array_events), keep.extra.columns = TRUE)
miso_events_granges <- makeGRangesFromDataFrame(as.data.frame(miso_bed), keep.extra.columns = TRUE)

##----------------------------------
## intersect files
##----------------------------------
hits <- findOverlaps(miso_events_granges, array_events_granges, minoverlap = 1, type = "any", select = "all")

eventKey <- data_frame(event_name = as.data.frame(miso_events_granges[queryHits(hits),])$name,
                       event_class = as.data.frame(miso_events_granges[queryHits(hits),])$score, 
                       gene_symbol = as.data.frame(array_events_granges[subjectHits(hits),])$name,
                       chr = as.data.frame(array_events_granges[subjectHits(hits),])$seqname,
                       start = as.data.frame(array_events_granges[subjectHits(hits),])$start,
                       end = as.data.frame(array_events_granges[subjectHits(hits),])$end,
                       Control_RTPCR = as.data.frame(array_events_granges[subjectHits(hits),])$Control_psi_mean,
                       DM1_RTPCR = as.data.frame(array_events_granges[subjectHits(hits),])$DM1_psi_mean,
                       array_event_class = as.data.frame(array_events_granges[subjectHits(hits),])$score)

eventKey <- unique(eventKey)
eventKey <- eventKey %>% filter(event_class %in% c("SE", "MXE"), array_event_class == "alt int ex")

write.table(eventKey, file = "~/Projects/DMseq/data/tibialis/array_intersect_miso_events.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
