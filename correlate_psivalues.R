##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param consolidated_summary_file 
##' @param pdata_file 
##' @param gff_file 
##' @param event_key_file 
##' @param grouping_var 
##' @param var1 
##' @param var2 
##' @param cor_method 
##' @param counts_filter 
##' @param inform_counts 
##' @param max_fraction_missing 
##' @return 
##' @author Adam Struck
calc_correlations <- function(consolidated_summary_file, pdata_file = NULL,
                              gff_file = NULL, event_key_file = NULL,
                              grouping_var = "isoforms", var1 = "miso_posterior_mean",
                              var2 = c("HG.QMT", "ADF.QMT", "Actual.Strength.6pt.Scale"),
                              cor_method = c("pearson", "pearson", "spearman"), exclude_samples = NULL,
                              counts_filter = TRUE, inform_counts = 20, max_fraction_missing = 0.25)
{
    require(data.table)
    require(dplyr)
    require(lazyeval)
    require(foreach)
    
    consolidated_summaries <- tbl_dt(fread(consolidated_summary_file))
    consolidated_summaries <- consolidated_summaries %>% filter(!sample %in% exclude_samples)
    n_samples <- consolidated_summaries$sample %>% unique() %>% length()

    if (counts_filter) {
        if (!"sum_informative_counts" %in% names(consolidated_summaries)) {
            source("./sum.informativeMisoCounts.R")
            consolidated_summaries[, sum_informative_counts := sum.InformativeMisoCounts(consolidated_summaries$counts)]
        }
        consolidated_summaries <- consolidated_summaries %>% filter(sum_informative_counts >= inform_counts)
    }

    if (!is.null(pdata_file)) {
        pdata <- fread(pdata_file)
        setnames(pdata, 1, "sample")
        pdata <- pdata %>% filter(!sample %in% exclude_samples, sample %in% unique(con_summaries$sample))
        ## combine consolidated data with pdata ##
        consolidated_summaries <- inner_join(consolidated_summaries, pdata, by = "sample")
    }
        
    ## Map MISO annotations to gene symbols ##
    if (!is.null(gff_file)) {
        gff <- fread(gff_file)
        gff <- gff %>% dplyr::filter(V3 == "gene")
        split <- strsplit(as.character(gff$V9), "[;=]")
        event.key <- data.table(do.call(rbind.data.frame, lapply(split, function(x) { 
            x[c(2,8,10,12)]
        })))
    }

    if (!is.null(event_key_file)) {
        event.key <- fread(event_key_file)
        setnames(event.key, c("event_name", "ensg_id", "gene_symbol", "desc"))
    }

    if(grouping_var == "isoforms") {
        named.events <- consolidated_summaries %>% select(event_name, isoforms) %>% distinct()
    }
    
    ## Correlate psi values with strength measurements ##
    con_data <- consolidated_summaries %>%
        group_by_(grouping_var) %>%
        mutate(n_missing = n_samples - n())

    if (length(cor_method) != length(var2))) {
        cor_method <- rep(cor_method[1], length(var2))
    }
    
    cor_res_list <- foreach(i = 1:length(var2), .combine = list, .multicombine = TRUE) %do%
    {
        summarize_formula <- list(interp(~cor.test(x, y, method = cor_method[i], use = "pairwise.complete.obs")$estimate, x = as.symbol(var1), y = as.symbol(var2[i])),
                                  interp(~cor.test(x, y, method = cor_method[i], use = "pairwise.complete.obs")$p.value, x = as.symbol(var1), y = as.symbol(var2[i])))
        res <- con_data %>%
            filter(n_missing / n_samples <= max_fraction_missing) %>%
            group_by_(grouping_var) %>%
            summarize_(.dots = summarize_formula)
        setnames(res, c(grouping_var, paste("cor", var1, var2[i], sep = "_"), paste("pval", var1, var2[i], sep = "_")))
        res
    }
    
    ## Merge results ##
    cor_data <- Reduce(function(...) merge(..., by = grouping_var, all = TRUE), cor_res_list)
    
    if(exists("named.events")) {
        cor_data <- merge(cor_data, named.events, by = grouping_var)
    }

    if(exists("event.key")) {
        cor_data <- merge(cor_data, event.key, by = "event_name", all.x = TRUE)
    }

    return(cor_data)
}
