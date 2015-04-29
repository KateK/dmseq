##----------------------------------
## required libraries
##----------------------------------
library(foreach)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

##----------------------------------
## event_type
##----------------------------------
chromosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                 "chr22", "chrM", "chrX", "chrY")


miso_events <- tbl_dt(fread("~/Projects/DMseq/data/event_annotations/nonUTRevents.multi.hg19.gff3"))

miso_events <- miso_events %>% dplyr::filter(V3 == "gene", V1 %in% chromosomes)

##----------------------------------
## transform miso annotations to bed
##----------------------------------
miso_events <- miso_events %>% select(V1, V4, V5, V9, V3, V7, V2)
setnames(miso_events, c("chr", "start", "end", "name", "score", "strand", "class"))

events <- gsub("ID=.+;", '', miso_events$name)
events <- gsub("Name=", "", events)

miso_bed <- do.call(rbind.data.frame, lapply(events, getEventCoords))
miso_bed <- tbl_df(miso_bed)

write.table(miso_bed, "~/Projects/DMseq/data/event_annotations/nonUTRevents.multi.hg19.bed", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

##----------------------------------
## function
##----------------------------------
getEventCoords <- function(event_str) {
    A3SS <- "^chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:|0-9]+[:][+-]$"
    A5SS <- "^chr[0-9XYM]+[:][0-9]+[:|0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][+-]$"
    SE <- "^chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]$"
    MXE <- "^chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]$"
    RI <- "^chr[0-9XYM]+[:][0-9]+[:|0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:|0-9]+[:][+-]$"
    if (str_count(event_str, A3SS) == 1) {
        event_type <- "A3SS"
        event_coord <- str_extract(event_str ,perl("(?<=@)chr[0-9XYM]+[:][0-9]+[:|0-9]+[:][+-]"))
        event_coord <- str_split(event_coord, "[:|]")
    } else {
        if (str_count(event_str, A5SS) == 1) {
            event_type <- "A5SS"
            event_coord <- str_extract(event_str ,perl("(?<!@)chr[0-9XYM]+[:][0-9]+[:|0-9]+[:][+-]"))
            event_coord <- str_split(event_coord, "[:|]")
        } else {
            if (str_count(event_str, SE) == 1) {
                event_type <- "SE"
                event_coord <- str_extract(event_str ,perl("(?<=@)chr[0-9XYM]{1,2}:[0-9]+:[0-9]+:[+-]"))
                event_coord <- str_split(event_coord, "[:|]")
            } else {
                if (str_count(event_str, MXE) == 1) {
                    event_type <- "MXE"
                    event_coord <- unlist(str_extract_all(event_str ,perl("(?<=@)chr[0-9XYM]{1,2}:[0-9]+:[0-9]+:[+-]")))
                    event_coord <- str_split(event_coord, "[:|]")
                    event_coord <- event_coord[1:2]
                } else {
                    if (str_count(event_str, RI) == 1) {
                        event_type <- "RI"
                        event_coord <- str_extract(event_str ,perl("(?<=@)chr[0-9XYM]{1,2}:[0-9]+[:|0-9]+:[+-]"))
                        event_coord <- str_split(event_coord, "[:|]")
                    } else {
                        return("Error: event didn't match any of the defined types")
                    }
                }
            }
        }
    }
    res <- do.call(rbind.data.frame, lapply(event_coord, function(event_coord)
        {
            data.frame(chr = event_coord[1],
                       start = event_coord[2],
                       end = event_coord[length(event_coord)-1],
                       name = event_str, score = event_type,
                       strand = event_coord[length(event_coord)])
        }))
    
    return(res)
}

