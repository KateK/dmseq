#!/usr/bin/env Rscript

##----------------------------------
## load analysis functions
##----------------------------------
setwd("~/Projects/DMseq/")

source("bin/summarizeMISOoutput_dev.R")

##----------------------------------
## sample group
##----------------------------------
sample_group <- "tibialis"
# sample_group <- "quadricep"
# sample_group <- "heart"
# sample_group <- "allControls"

##----------------------------------
## base directories
##----------------------------------
base_dir <- paste("data/", sample_group, "/", sep = "")
# miso_dir <- paste("~/MISO/single_end_mode/summaries/", sample_group, "/", sep = "")
miso_dir <- base_dir
output_dir <- paste("results/", sample_group, "/", sep = "")

##----------------------------------
## phenotypic data about samples
##----------------------------------
# pdata_file <- paste("data/", paste("DM", sample_group, "pdata.txt", sep = "_"), sep = "")
pdata_file <- "data/DM_sample_pdata.txt"

##----------------------------------
## excluded samples
##----------------------------------
# library(data.table)
# pdata <- fread(pdata_file)
## pdata

# excluded_samples <- c("R153")
# excluded_samples <- c("BCM270")
# excluded_samples <-c("MN090909_card")
# excluded_samples <-c("BCM272")
# excluded_samples <- NULL

##----------------------------------
## all event types
##----------------------------------
event_types <- c("nonUTRevents.multi", "ALE", "AFE", "polyA")
## event_types <- c("polyA")
compare_by <- c("Diagnosis")

for (i in 1:length(event_types)) {
    event_type <- event_types[i]
    message("Working on: ", event_type, Sys.time())
    consolidated_summary_file <- paste(miso_dir, paste(sample_group, event_type, "consolidatedSummaries.txt", sep = "_"), sep = "")
    consolidated_bayesfactor_file <- paste(miso_dir, paste(sample_group, event_type, "consolidatedBayesFactors.txt", sep = "_"), sep = "")
    output_file <- paste(output_dir, paste(sample_group, event_type, "results.txt", sep = "_"), sep = "")
    if (event_type == "nonUTRevents.multi") {
        # event_key_file <- "~/alt.splicing.annotations/hg19_Aug2014/nonUTRevents.multi.hg19.to.ensGene.txt"
        event_key_file <- "data/event_annotations/nonUTRevents.multi.hg19.to.ensGene.txt"
    } else {
        event_key_file <- NULL
    }
    
    res_data <- summarizeMISOoutput(consolidated_summary_file, consolidated_bayesfactor_file, pdata_file,
                                    main_grouping_var = "isoforms", aggregate_stats_by_vars = compare_by,
                                    aggregate_stats_for_vars = c("miso_posterior_mean"),
                                    bf_wide_format = FALSE, counts_filter = TRUE, inform_counts = 20,
                                    event_key_file = event_key_file, exclude_samples = excluded_samples)
    
    write.table(res_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    message("Finished: ", event_type, Sys.time())
}

## ##----------------------------------
## ## nonUTRevents.multi
## ##----------------------------------
## consolidated_summary_file <- paste(miso_dir, paste(sample_group, "nonUTRevents.multi_consolidatedSummaries.txt", sep = "_"), sep = "")
## consolidated_bayesfactor_file <- paste(miso_dir, paste(sample_group, "nonUTRevents.multi_consolidatedBayesFactors.txt", sep = "_"), sep = "")
## event_key_file <- "~/alt.splicing.annotations/hg19_Aug2014/nonUTRevents.multi.hg19.to.ensGene.txt"
## output_file <- paste(output_dir, paste(sample_group, "nonUTRevents.multi_results.txt", sep = "_"), sep = "")

## res_data <- summarizeMISOoutput(consolidated_summary_file, consolidated_bayesfactor_file, pdata_file = pdata_file,
##                                 main_grouping_var = "isoforms", bf_wide_format = FALSE, compare_strength_groups = FALSE,
##                                 counts_filter = TRUE, event_key_file = event_key_file, exclude_samples = excluded_samples)

## write.table(res_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## ##----------------------------------
## ## ALE
## ##----------------------------------
## consolidated_summary_file <- paste(miso_dir, paste(sample_group, "ALE_consolidatedSummaries.txt", sep = "_"), sep = "")
## consolidated_bayesfactor_file <- paste(miso_dir, paste(sample_group, "ALE_consolidatedBayesFactors.txt", sep = "_"), sep = "")
## event_key_file <- NULL
## output_file <- paste(output_dir, paste(sample_group, "ALE_results.txt", sep = "_"), sep = "")

## res_data <- summarizeMISOoutput(consolidated_summary_file, consolidated_bayesfactor_file, pdata_file = pdata_file,
##                                 main_grouping_var = "isoforms", bf_wide_format = FALSE, compare_strength_groups = FALSE,
##                                 counts_filter = TRUE, event_key_file = event_key_file, exclude_samples = excluded_samples)

## write.table(res_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

## ##----------------------------------
## ## AFE
## ##----------------------------------
## consolidated_summary_file <- paste(miso_dir, paste(sample_group, "AFE_consolidatedSummaries.txt", sep = "_"), sep = "")
## consolidated_bayesfactor_file <- paste(miso_dir, paste(sample_group, "AFE_consolidatedBayesFactors.txt", sep = "_"), sep = "")
## event_key_file <- NULL
## output_file <- paste(output_dir, paste(sample_group, "AFE_results.txt", sep = "_"), sep = "")

## res_data <- summarizeMISOoutput(consolidated_summary_file, consolidated_bayesfactor_file, pdata_file = pdata_file,
##                                 main_grouping_var = "isoforms", bf_wide_format = FALSE, compare_strength_groups = FALSE,
##                                 counts_filter = TRUE, event_key_file = event_key_file, exclude_samples = excluded_samples)

## write.table(res_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
## CLCN1_multi
##----------------------------------
## consolidated_summary_file <- paste(miso_dir, paste(sample_group, "CLCN1_multi_consolidatedSummaries.txt", sep = "_"), sep = "")
## consolidated_bayesfactor_file <- paste(miso_dir, paste(sample_group, "CLCN1_multi_consolidatedBayesFactors.txt", sep = "_"), sep = "")
## event_key_file <- NULL
## output_file <- paste(output_dir, paste(sample_group, "CLCN1_multi_results.txt", sep = "_"), sep = "")

## res_data <- summarizeMISOoutput(consolidated_summary_file, consolidated_bayesfactor_file, pdata_file = pdata_file,
##                                 main_grouping_var = "isoforms", bf_wide_format = FALSE, compare_strength_groups = TRUE,
##                                 counts_filter = TRUE, event_key_file = event_key_file, exclude_samples = excluded_samples)

## write.table(res_data, file = output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
