#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

summary_dir  <- as.character(args[1])
file_type    <- as.character(args[2])
header       <- as.logical(args[3])
out          <- args[4]

message("Output file: ", out, "\n")

analyzeMicroExonEvents <- function(summary_dir, file_type = "sam", header = TRUE, output_file) {

    require(rbamtools)
    require(foreach)
    require(dplyr)
    require(tidyr)
    require(data.table)
    
    files <- list.files(path = summary_dir, pattern = file_type, full.names = TRUE)
    if (length(files) == 0) stop("ERROR: incorrect summary directory or file extension provided")

    output_table <- foreach(i = 1:length(files), .combine = rbind, .multicombine = TRUE) %do%
    {
        file <- files[i]
        message("Working on file: ", i, "/", length(files), " ", Sys.time())
        filename       <- basename(file)
        filename_base  <- unlist(strsplit(filename, "[.]"))[1]

        reader <- bamReader(file, idx=TRUE)
        counts <- bamCountAll(reader, verbose = FALSE)
        bamClose(reader)
        
        counts$event_name <- rownames(counts)
        rownames(counts) <- NULL
        data.frame(sample = filename_base, event_name = counts$event_name, nAligns = counts$M, LN = counts$LN)
    }

    output_table <- tbl_df(output_table)
    event_info <- splitMicroExonID(output_table$event_name)

    mic_info <- tbl_df(cbind(output_table, event_info))
    
    mic_counts <- spread(mic_info, value = nAligns, key = event_class) %>% select(sample, event_name, gene_id)
    
}


splitMicroExonID <- function(x) {
    require(stringr)
    event_info <- str_split(x, "[-.]")
    do.call(rbind.data.frame, lapply(event_info, function(event) {
        data.frame(gene_id = event[1], event_class = event[length(event) - 1])
    }))
}
