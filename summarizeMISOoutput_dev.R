##' Summarize MISO analysis output
##'
##' Calculate summary statistics [mean & sd] for each event for defined groups
##' 
##' @param consolidated_summary_file consolidated sample summaries from MISO 
##' @param consolidated_bayesfactor_file [optional] consolidated sample comparisons from MISO
##' @param pdata_file metadata about samples
##' @param main_grouping_var if MISO was run using multi-isoform annotations use "isoforms", otherwise use "event_name"
##' @param aggregate_stats_by_vars [default: "Diagnosis"] aggregate statistics using these sample grouping variables
##' @param aggregate_stats_for_vars [default: "miso_posterior_mean"] aggregate statistics for these sample variables
##' @param counts_filter [TRUE/FALSE] filter the data based on the number of informative reads over each event
##' @param inform_counts [default: 20] filtering criteria - minimum number of informative inclusion/exclusion reads required for an event
##' @param exclude_samples [optional] which samples (if any) should be excluded from the analysis?
##' @param event_key_file [optional] used to map 'main_grouping_var' to gene_symbol
##' @param bf_wide_format [TRUE/FALSE] is consolidated_bayesfactor_file in wide format?
##' @return data.frame / data.table 
##' @author Adam Struck
summarizeMISOoutput <- function(consolidated_summary_file, consolidated_bayesfactor_file = NULL,
                                pdata_file, main_grouping_var = "isoforms", aggregate_stats_by_vars = c("Diagnosis"),
                                aggregate_stats_for_vars = c("miso_posterior_mean"), 
                                counts_filter = TRUE, inform_counts = 20,
                                exclude_samples = NULL, event_key_file = NULL)
{
    require(data.table)
    require(dplyr)
    require(tidyr)
    require(foreach)
    require(stringr)
    
    ### Load files into memory ###
    con_summaries <- tbl_dt(fread(consolidated_summary_file))
    con_summaries[, ci_width := ci_high - ci_low]
    con_summaries <- con_summaries %>% filter(!sample %in% exclude_samples)
    
    if (counts_filter) {
        if (!"sum_informative_counts" %in% names(con_summaries)) {
            source("./sum.informativeMisoCounts.R")
            con_summaries[, sum_informative_counts := sum.InformativeMisoCounts(con_summaries$counts)]
        }
        con_summaries <- con_summaries %>% filter(sum_informative_counts >= inform_counts)
    }

    pdata <- fread(pdata_file, header = TRUE) 
    setnames(pdata, 1, "sample")
    pdata <- pdata %>% filter(!sample %in% exclude_samples, sample %in% unique(con_summaries$sample))
    
    if (!is.null(event_key_file)) {
        event.key <- fread(event_key_file)
        setnames(event.key, c("event_name", "ensg_id", "gene_symbol", "desc"))
    }
    
    if(main_grouping_var == "isoforms") {
        named.events <- con_summaries %>% select(event_name, isoforms) %>% distinct()
    }
    
    ## combine consolidated data with pdata
    consolidated_data <- inner_join(con_summaries, pdata, by = "sample")

    ## summarize by grouping vars
    all_summaries <- foreach(i = 1:length(aggregate_stats_for_vars), .combine = c, .multicombine = TRUE) %do% {
        aggregate_stats_for_var <- aggregate_stats_for_vars[i]
        if (i == 1) { return_n_logical <- TRUE } else { return_n_logical <- FALSE }
        foreach(j = 1:length(aggregate_stats_by_vars), .combine = c, .multicombine = TRUE) %do% {
            aggregate_stats_by_var <- aggregate_stats_by_vars[j]
            summarize_field(consolidated_data, field = aggregate_stats_for_var,
                            grouping_vars = c(main_grouping_var, aggregate_stats_by_var),
                            newname = aggregate_stats_for_var, return_n = return_n_logical)
        }
    }

    ## combine all group summaries into data table
    results_summary <- Reduce(function(...) merge(..., by = main_grouping_var, all = TRUE), all_summaries)

    ## summarize comparisons
    if (!is.null(consolidated_bayesfactor_file)) {
        con_bf <- tbl_dt(fread(consolidated_bayesfactor_file))
        con_bf <- con_bf %>% filter(!sample1 %in% excluded_samples,!sample2 %in% excluded_samples)
        if (counts_filter) {
            if (bf_wide_format) {
                return("WARNING: Data was not filtered. Use a consolidated_bayesfactor_file in long format.")
            } else {
                x <- paste(con_summaries$sample, con_summaries[[main_grouping_var]], sep ="_")
                y <- paste(con_bf$sample1, con_bf[[main_grouping_var]], sep ="_")
                z <- paste(con_bf$sample2, con_bf[[main_grouping_var]], sep ="_")
                con_bf <- con_bf[y %in% x & z %in% x,]
                rm(list = c("x", "y", "z"))
            }
        }
        
        all_n_sig <- foreach(i = 1:length(aggregate_stats_by_vars), .combine = list, .multicombine = TRUE) %do% {
            aggregate_stats_by_var <- aggregate_stats_by_vars[i]
            group_comps <- get_group_comps(aggregate_stats_by_var, pdata)
            n_sig <- summarizeBF(group_comps,
                                 grouping_var = aggregate_stats_by_var,
                                 main_grouping_var, con_bf, pdata)
        }
        
        if (!is.data.table(all_n_sig)) {
            ## names(all_n_sig) <- aggregate_stats_by_vars
            n_sig <- Reduce(function(...) merge(..., by = main_grouping_var, all = TRUE), all_n_sig)
        }
        
        results_summary <- merge(results_summary, n_sig, by = main_grouping_var, all = TRUE)    
    }

    if (exists("named.events")) {
        results_summary <- left_join(results_summary, named.events, by = "isoforms")
    }
    
    if (exists("event.key")) {
        results_summary <- left_join(results_summary, event.key, by = "event_name")
    }
    
    results_summary <- results_summary[, c(sort(names(results_summary))), with = FALSE]
    return(results_summary)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title summarize_field
##' @param inputDT data.table object of consolidated summaries
##' @param field variable for which summary stats will be calculated
##' @param grouping_vars What variables should observations be grouped by prior to aggregation?
##' @param newname new name for "field"
##' @param return_n [TRUE/FALSE] return number of non-NA observations for "field"
##' @return data.table
##' @author Adam Struck
summarize_field <- function(inputDT, field, grouping_vars, newname, return_n = TRUE) {
    if (is.character(grouping_vars) & length(grouping_vars) != 2) {
        stop("ERROR: grouping_vars must be a character vector of length 2")
    }
    if (field == "miso_posterior_mean") {
        newname <- "psi"
    }
    ## setup aggreagation functions
    field_summary <- list(mean = ~mean(x, na.rm=TRUE), sd = ~sd(x, na.rm=TRUE), n = ~sum(!is.na(x)))
    field_summary <- lapply(field_summary, lazyeval::interp, x = as.symbol(field)) 
    ## aggregate stats
    int_res <- inputDT %>%
        group_by_(grouping_vars[1], grouping_vars[2]) %>%
        summarize_(.dots = field_summary)
    ## split into wide format
    stats <- c("mean", "sd")
    for (i in 1:length(stats)) {
        stat <- stats[i]
        assign(paste(newname, stat, "by", grouping_vars[2], sep = "_"),
               int_res %>% ungroup %>% select_(grouping_vars[1], grouping_vars[2], stat) %>% spread_(grouping_vars[2], stat))
        setnames(get(paste(newname, stat, "by", grouping_vars[2], sep = "_")),
                 c(grouping_vars[1], paste(sort(unique(inputDT[[grouping_vars[2]]])), newname, stat, sep = "_")))
    }
    if (return_n) {
        assign(paste("n", "by", grouping_vars[2], sep = "_"),
               int_res %>% ungroup %>% select_(grouping_vars[1], grouping_vars[2], "n") %>% spread_(grouping_vars[2], "n"))
        setnames(get(paste("n", "by", grouping_vars[2], sep = "_")),
                 c(grouping_vars[1], paste(sort(unique(inputDT[[grouping_vars[2]]])), "n", sep = "_")))        
        res <- list(get(paste(newname, "mean", "by", grouping_vars[2], sep = "_")),
                    get(paste(newname, "sd", "by", grouping_vars[2], sep = "_")),
                    get(paste("n", "by", grouping_vars[2], sep = "_")))
        names(res) <- c(paste(newname, "mean", "by", grouping_vars[2], sep = "_"),
                        paste(newname, "sd", "by", grouping_vars[2], sep = "_"),
                        paste("n", "by", grouping_vars[2], sep = "_"))
    } else {
        res <- list(get(paste(newname, "mean", "by", grouping_vars[2], sep = "_")),
                    get(paste(newname, "sd", "by", grouping_vars[2], sep = "_")))
        names(res) <- c(paste(newname, "mean", "by", grouping_vars[2], sep = "_"),
                        paste(newname, "sd", "by", grouping_vars[2], sep = "_"))
    }
    return(res)
}

##' Get all possible sample comparisons for set groups
##'
##' .. content for \details{} ..
##' @title get_group_comps
##' @param group_var variable to group samples by
##' @param pdata data.table object containing sample metadata
##' @return list
##' @author Adam Struck
get_group_comps <- function(group_var, pdata) {
    groups <- unique(pdata[[group_var]])
    group_comps <- NULL
    res <- foreach(i = 1:(length(groups)-1), .combine = c, .multicombine = TRUE) %do%
    {
        group1 <- groups[i]
        group1_samples <-  pdata %>% filter_(.dots = list(lazyeval::interp(~x == group1, x = as.symbol(group_var))))
        group1_samples <-  group1_samples$sample
        foreach(j = (i+1):length(groups), .combine = c, .multicombine = TRUE) %do%
        {
            group2 <- groups[j]
            group_comps <- c(group_comps, paste(group1, group2, sep = "_vs_"))
            group2_samples <-  pdata %>% filter_(.dots = list(lazyeval::interp(~x == group2, x = as.symbol(group_var))))
            group2_samples <- group2_samples$sample
            out <- NULL
            res <- foreach(k = 1:length(group2_samples), .combine = c, .multicombine = TRUE) %do%
            {
                c(paste(group1_samples, group2_samples[k], sep = "_vs_"), paste(group2_samples[k], group1_samples, sep = "_vs_"))
            }
            list(res)
        }
    }
    names(res) <- group_comps
    return(res)
}

##' Summarize MISO comparison data
##'
##' .. content for \details{} ..
##' @title summarizeBF
##' @param group_comps list of sample comparisons
##' @param grouping_var variable that seperates the samples are we comparing into groups
##' @param main_grouping_var how events are grouped; isoform or event (some events have more than 2 isoforms)
##' @param con_bf data.table with columns: sample1, sample2, diff, bayes_factor
##' @param pdata data.table with sample information for grouping
##' @param threshold list of thresholds to use when summarizing group-level comparisons 
##' @return data.table
##' @author Adam Struck
summarizeBF <- function(group_comps, grouping_var, main_grouping_var, con_bf, pdata, threshold = list(deltapsi = 0.05, bf = 5, proportion = 0.5)) {
    if (!"comp" %in% names(con_bf)) {
        con_bf[, comp := paste(sample1, "vs", sample2, sep = "_")]
    }
    n_sig <- foreach(i = 1:length(group_comps), .combine = list, .multicombine = TRUE) %do%
    {           
        group1 <- unlist(str_split(names(group_comps)[i], "_vs_"))[1]
        group1_samples <- pdata[pdata[[grouping_var]] == group1,]$sample
        group2 <- unlist(str_split(names(group_comps)[i], "_vs_"))[2]
        group2_samples <- pdata[pdata[[grouping_var]] == group2,]$sample

        n_comps <- con_bf %>%
            select(sample1, sample2, comp) %>%
            filter(comp %in% group_comps[[i]]) %>% 
            distinct %>%
            group_by(sample1) %>%
            summarize(n_comps = n())
        n_comp_name <- paste(names(group_comps)[i], "n_comps", sep = "_")
        setnames(n_comps, c("sample1", n_comp_name))
        
        tmp_n_sig <- con_bf %>%
            filter(sample1 %in% group1_samples & sample2 %in% group2_samples | sample2 %in% group1_samples & sample1 %in% group2_samples,
                   comp %in% group_comps[[i]],
                   abs(diff) >= threshold$deltapsi,
                   bayes_factor >= threshold$bf) %>%
            group_by_(main_grouping_var, "sample1") %>%            
            summarize(n_sig = n())
        n_sig_name <- paste(names(group_comps)[i], "n_sig", sep = "_")
        setnames(tmp_n_sig, c(main_grouping_var, "sample1", n_sig_name))

        tmp_n_sig <- left_join(n_comps, tmp_n_sig, by = "sample1")

        filter_criteria <- interp(~x / y > threshold$proportion,
                                  x = as.symbol(n_sig_name),
                                  y = as.symbol(n_comp_name))
        n_sig_res <- tmp_n_sig %>%
            filter_(.dots = list(filter_criteria)) %>%
            group_by_(main_grouping_var) %>%
            summarize(n_sig = n())
        setnames(n_sig_res, c(main_grouping_var, n_sig_name))
    }
    n_sig <- Reduce(function(...) merge(..., by = main_grouping_var, all = TRUE), n_sig)
    return(tbl_dt(n_sig))
}
