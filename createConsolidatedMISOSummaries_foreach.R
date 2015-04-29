#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

summary_dir  <- as.character(args[1])
file_type    <- as.character(args[2])
header       <- as.logical(args[3])
out          <- args[4]

message("Output file: ", out, "\n")

create_consolidated_miso_summary <- function(summary_dir, file_type = "miso_summary", header = TRUE, output_file) {

    # require(doMC)
    # registerDoMC(cores = detectCores())
    require(foreach)
    require(data.table)
    
    files <- list.files(path = summary_dir, pattern = file_type, full.names = TRUE)
    if (length(files) == 0) stop("ERROR: incorrect summary directory or file extension provided")

    comps <- vector()
    output_table <- foreach(i = 1:length(files), .packages = c("data.table"), .combine = rbind, .multicombine = TRUE) %do%
    {
        file <- files[i]
        message("Working on file: ", i, "/", length(files), " ", Sys.time())
        filename       <- basename(file)
        filename_ext   <- unlist(strsplit(filename, "[.]"))[2]
        filename_base  <- unlist(strsplit(filename, "[.]"))[1]
        
        if ( filename_ext == "miso_summary" ) {
            data_f <- as.data.table(read.table(file = file, sep = "\t", header = header, stringsAsFactors = FALSE))
            attributes <- c("event_name", "ci_low", "miso_posterior_mean", "ci_high", "counts", "mRNA_starts", "mRNA_ends")
            data <- data_f[, attributes, with = FALSE]
            ## check if multiiso 
            if (dim(data)[1] > 0) {
                if (!is.numeric(data$ci_low)) {
                    data <- data.table(data, data_f[, "isoforms", with = FALSE])
                    ## split values in columns on "," into new rows
                    data <- split_csv_cols(data, filename_ext)
                }
                data <- data.table(sample = filename_base, data)
            }
        } else {
            if ( filename_ext == "miso_bf" ) {
                sample1_base   <- unlist(strsplit(filename_base, "_vs_"))[1]
                sample2_base <- unlist(strsplit(filename_base, "_vs_"))[2]
                if (sample1_base != sample2_base) {
                    if (!filename_base %in% comps) {
                        comps <- c(comps, paste(sample1_base, "_vs_", sample2_base, sep = ""), paste(sample2_base, "_vs_", sample1_base, sep = ""))
                        attributes <- c("event_name", "diff", "bayes_factor")
                        data_f <- as.data.table(read.table(file = file, sep = "\t", header = header, stringsAsFactors = FALSE))
                        data <- data_f[, attributes, with = FALSE]
                        ## check if multi iso
                        if (dim(data)[1] > 0) {
                            if (!is.numeric(data$bayes_factor)) {
                                data <- data.table(data, data_f[, "isoforms", with = FALSE])
                                data <- split_csv_cols(data, filename_ext)
                            }
                            data <- data.table(sample1 = sample1_base,
                                               sample2 = sample2_base,
                                               comp = paste(sample1_base, "_vs_", sample2_base, sep=""),
                                               data)
                        }
                    }
                }
            }
        }
        message("Finished processing ", file, " ", Sys.time())
        if(dim(data)[1] > 0) {
            data
        }
    }
    write_results(output_table, output_file)
}

### subroutines ###
split_counts_field <- function(counts_field) {
    categories <- unlist(strsplit(counts_field, "[:][0-9]+[,]?"))
    counts <- unlist(
        lapply(strsplit(counts_field, "[(][01,]+[)][:]"), function(x) {
            as.numeric(sub(",", "", x))
        }))
    counts <- setNames(counts[-1], categories)
    result <- list(reads_discarded = counts["(0,0)"], reads_first_iso = counts["(1,0)"],
                   reads_second_iso = counts["(0,1)"], reads_both_iso = counts["(1,1)"])
    return(result)
}

sum.InformativeMisoCounts <- function(counts_string) {
    ## remove uninformative counts
    counts <- gsub("\\([0,]+\\):[0-9]+[,]?", "", counts_string, perl = TRUE)
    counts <- gsub("[,]?\\([1,]+\\):[0-9]+", "", counts, perl = TRUE)
    ## remove counts class identifier and sum counts
    ## of informative reads
    inform_counts <- gsub("\\([01,]+\\):", "", counts, perl = TRUE)
    sapply(str_split(inform_counts, ","), function(x) {
        sum(as.numeric(x))
    })
}

splitter <- function(dt, x, ul = TRUE) {
    if (ul) {
        tmp <- strsplit(dt[[x]], ",")
        unlist(sapply(tmp, function(x) {
            if (length(x) == 2) {
                paste(x, collapse = ",")
            } else {
                x
            }
        }))
    } else {
        tmp <- strsplit(dt[[x]], ",")
        sapply(tmp, function(x) {
            if (length(x) == 2) {
                paste(x, collapse = ",")
            } else {
                x
            }
        })
    }
}

split_csv_cols <- function(dt, file_type) {
    if (file_type == "miso_summary") {
        res <- data.table(event_name = rep(dt$event_name, sapply(splitter(dt, "ci_low", ul = FALSE), length)),
                          ci_low = splitter(dt, "ci_low"),
                          miso_posterior_mean = splitter(dt, "miso_posterior_mean"),
                          ci_high = splitter(dt, "ci_high"),
                          sum_informative_counts = rep(sum.InformativeMisoCounts(dt$counts), sapply(splitter(dt, "ci_low", ul = FALSE), length)),
                          assigned_counts = splitter(dt, "assigned_counts"),
                          chr = rep(dt$chrom, sapply(splitter(dt, "ci_low", ul = FALSE), length)),
                          mRNA_starts = splitter(dt, "mRNA_starts"),
                          mRNA_ends = splitter(dt, "mRNA_ends"),
                          isoforms = splitter(dt, "isoforms"))
    } else {
        if (file_type == "miso_bf") {
            res <- data.table(event_name = rep(dt$event_name, sapply(splitter(dt, "bayes_factor", ul = FALSE), length)),
                              isoforms = splitter(dt, "isoforms"),
                              diff = splitter(dt, "diff"),
                              bayes_factor = splitter(dt, "bayes_factor"))
        }
    }
    return(res)
}

write_results <- function(results.table, output_file) {
    write.table(results.table,
                col.names = TRUE,
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                file = output_file)
}

create_consolidated_miso_summary(summary_dir = summary_dir,
                                     file_type = file_type,
                                     header = header,
                                     output_file = out)
