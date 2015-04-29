##' Summarize MISO analysis output
##'
##' Calculate summary statistics for each event grouped by 'Diagnosis' and optionally, 'Strength_Grouping'
##' 
##' @param consolidated_summary_file 
##' @param pdata_file metadata about samples
##' @param consolidated_bayesfactor_file [optional] used to summarizes the number of samples that were, on average, different from controls
##' @param bf_wide_format [TRUE/FALSE] is consolidated_bayesfactor_file in wide format?
##' @param main_grouping_var if MISO was run using multi-isoform annotations use "isoforms", otherwise use "event_name"
##' @param compare_strength_groups [TRUE/FALSE] summarize grouping by 'Strength_Grouping' data in pdata file
##' @param counts_filter [TRUE/FALSE] filter the data based on the number of informative reads over each event
##' @param inform_counts [default: 20] filtering criteria - minimum number of informative inclusion/exclusion reads required for an event
##' @param exclude_samples [optional] which samples (if any) should be excluded from the analysis?
##' @param event_key_file [optional] used to map 'main_grouping_var' to gene_symbol
##' @return data.frame / data.table 
##' @author Adam Struck
summarizeMISOoutput <- function(consolidated_summary_file, pdata_file,
                                consolidated_bayesfactor_file = NULL, bf_wide_format = FALSE,
                                main_grouping_var = "isoforms", compare_strength_groups = TRUE,
                                counts_filter = TRUE, inform_counts = 20, exclude_samples = NULL,
                                event_key_file = NULL) {

    require(data.table)
    require(dplyr)
    require(reshape2)
    require(foreach)
    require(stringr)
    
    ### Load files into memory ###
    con_summaries <- tbl_dt(fread(consolidated_summary_file))
    con_summaries[, ci_width := ci_high - ci_low]
    con_summaries <- con_summaries %>% filter(!sample %in% exclude_samples)
    
    if (counts_filter) {
        con_summaries[, informativeCounts := sum.InformativeMisoCounts(con_summaries$counts)]
        con_summaries <- con_summaries[informativeCounts >= inform_counts,]        
    }

    pdata <- read.table(pdata_file, header = TRUE, row.names = NULL, sep = "\t") 
    names(pdata)[1] <- "sample"
    pdata <- pdata %>% filter(!sample %in% exclude_samples)
    
    if (!is.null(event_key_file)) {
        event.key <- fread(event_key_file)
        setnames(event.key, c("event_name", "ensg_id", "gene_symbol", "desc"))
    }
    
    if(main_grouping_var == "isoforms") {
        named.events <- con_summaries %>% select(event_name, isoforms) %>% distinct()
    }
    
    ## combine consolidated data with pdata
    consolidated_data <- tbl_dt(merge(con_summaries, pdata, by = "sample", all = FALSE))

    psi_diagnosis_summary <- summarize_field(consolidated_data, field = "miso_posterior_mean",
                                             grouping_vars = c(main_grouping_var, "Diagnosis"),
                                             newname = "psi", return_n = TRUE)

    ci_diagnosis_summary <- summarize_field(consolidated_data, field = "ci_width",
                                            grouping_vars = c(main_grouping_var, "Diagnosis"),
                                            newname = "ci_width", return_n = TRUE)

    if (compare_strength_groups) {
        psi_strength_summary <- summarize_field(consolidated_data, field = "miso_posterior_mean",
                                                grouping_vars = c(main_grouping_var, "Strength_Grouping"),
                                                newname = "psi", return_n = TRUE)
    
        results_summary <- Reduce(function(...) merge(..., by = main_grouping_var, all = TRUE),
                                  list(psi_diagnosis_summary$psi_mean_by_Diagnosis,
                                       psi_diagnosis_summary$psi_sd_by_Diagnosis,
                                       psi_diagnosis_summary$n_by_Diagnosis,
                                       ci_diagnosis_summary$ci_width_mean_by_Diagnosis,
                                       ci_diagnosis_summary$ci_width_sd_by_Diagnosis,
                                       psi_strength_summary$psi_mean_by_Strength_Grouping[,c(1,3,4,5), with = F],
                                       psi_strength_summary$psi_sd_by_Strength_Grouping[, c(1,3,4,5), with = F],
                                       psi_strength_summary$n_by_Strength_Grouping[, c(1,3,4,5), with=F]))
    } else {
        results_summary <- Reduce(function(...) merge(..., by = main_grouping_var, all = TRUE),
                                  list(psi_diagnosis_summary$psi_mean_by_Diagnosis,
                                       psi_diagnosis_summary$psi_sd_by_Diagnosis,
                                       psi_diagnosis_summary$n_by_Diagnosis,
                                       ci_diagnosis_summary$ci_width_mean_by_Diagnosis,
                                       ci_diagnosis_summary$ci_width_sd_by_Diagnosis))
    }    

    if (!is.null(consolidated_bayesfactor_file)) {
        con_bf <- tbl_dt(fread(consolidated_bayesfactor_file))
        if (counts_filter) {
            if (bf_wide_format) {
                return("WARNING: Data was not filtered. Use a consolidated_bayesfactor_file in long format.")
            } else {
                x <- paste(con_summaries$sample, con_summaries[["isoforms"]], sep ="_")
                y <- paste(con_bf$sample1, con_bf[["isoforms"]], sep ="_")
                z <- paste(con_bf$sample2, con_bf[["isoforms"]], sep ="_")
                con_bf <- con_bf[y %in% x & z %in% x,]
                rm(list = c("x", "y", "z"))
            }
        }
        ## oragnize relevant comparisons into groups
        diagnosis_group_comps <- get_group_comps("Diagnosis", pdata)
        
        ## summarize sample comparisons by group
        n_sig <- summarizeBF(group_comps = diagnosis_group_comps,
                             grouping_var = "Diagnosis",
                             main_grouping_var, con_bf, pdata, bf_wide_format)

        if (compare_strength_groups) {
            strength_group_comps <- get_group_comps("Strength_Grouping", pdata)
            strength_group_comps$Control <- NULL
            
            n_sig_byStrength <- summarizeBF(group_comps = strength_group_comps,
                                            grouping_var = "Strength_Grouping",
                                            main_grouping_var, con_bf, pdata, bf_wide_format)

            n_sig <- merge(n_sig, n_sig_byStrength, by = main_grouping_var, all = TRUE)
        }
        results_summary <- merge(results_summary, n_sig, by = main_grouping_var, all = TRUE)
    }

    if (exists("named.events")) {
        results_summary <- merge(results_summary, named.events, by = "isoforms", all.x = TRUE)
    }
    
    if (exists("event.key")) {
        results_summary <- merge(results_summary, event.key, by = "event_name", all.x = TRUE)
    }
    
    results_summary[, delta_psi_mean := DM1_psi_mean - Control_psi_mean]
    results_summary <- results_summary[, c(sort(names(results_summary))), with = FALSE]
    results_summary <- results_summary[order(results_summary$DM1_n_sig, decreasing = TRUE),]
    return(results_summary)
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

summarize_field <- function(x, field, grouping_vars, newname, return_n = TRUE) {
    if (is.character(grouping_vars) & length(grouping_vars) != 2) {
        stop("ERROR: grouping_vars must be a character vector of length 2")
    }
    field_summary <- list(mean = ~mean(x, na.rm=TRUE), sd = ~sd(x, na.rm=TRUE), n = ~sum(!is.na(x)))
    field_summary <- lapply(field_summary, lazyeval::interp, x = as.symbol(field)) 
    tmp <- x %>% group_by_(grouping_vars[1], grouping_vars[2]) %>% summarize_(.dots = field_summary)
    ## field_summary <- list(mean = ~mean(get(field), na.rm=TRUE), sd = ~sd(get(field), na.rm=TRUE), n = ~sum(!is.na(get(field))))
    ## tmp <- x %>% group_by_(grouping_vars[1], grouping_vars[2]) %>% summarize(mean=mean(get(field), na.rm=TRUE), sd=sd(get(field), na.rm=TRUE), n=sum(!is.na(get(field))))

    assign(paste(newname, "mean", "by", grouping_vars[2], sep="_"), tbl_dt(dcast(tmp, get(grouping_vars[1]) ~ get(grouping_vars[2]), value.var = 'mean')))
    setnames(get(paste(newname, "mean", "by", grouping_vars[2], sep="_")), c(grouping_vars[1], paste(sort(unique(x[[grouping_vars[2]]])), newname, "mean", sep = "_")))

    assign(paste(newname, "sd", "by", grouping_vars[2], sep="_"), tbl_dt(dcast(tmp, get(grouping_vars[1]) ~ get(grouping_vars[2]), value.var = 'sd')))
    setnames(get(paste(newname, "sd", "by", grouping_vars[2], sep="_")), c(grouping_vars[1], paste(sort(unique(x[[grouping_vars[2]]])), newname, "sd", sep = "_")))

    assign(paste("n", "by", grouping_vars[2], sep="_"), tbl_dt(dcast(tmp, get(grouping_vars[1]) ~ get(grouping_vars[2]), value.var = 'n')))
    setnames(get(paste("n", "by", grouping_vars[2], sep="_")), c(grouping_vars[1], paste(sort(unique(x[[grouping_vars[2]]])), "n", sep = "_")))
    
    if (return_n) {
        res <- list(get(paste(newname, "mean", "by", grouping_vars[2], sep="_")), get(paste(newname, "sd", "by", grouping_vars[2], sep="_")), get(paste("n", "by", grouping_vars[2], sep="_")))
        names(res) <- c(paste(newname, "mean", "by", grouping_vars[2], sep="_"), paste(newname, "sd", "by", grouping_vars[2], sep="_"), paste("n", "by", grouping_vars[2], sep="_"))
    } else {
        res <- list(get(paste(newname, "mean", "by", grouping_vars[2], sep="_")), get(paste(newname, "sd", "by", grouping_vars[2], sep="_")))
        names(res) <- c(paste(newname, "mean", "by", grouping_vars[2], sep="_"), paste(newname, "sd", "by", grouping_vars[2], sep="_"))
    }
    return(res)
}    

## compare all levels of [group_var] to "Control" level
get_group_comps <- function(group_var, pdata) {
    groups <- unique(pdata[, group_var])
    group_comps <- NULL
    res <- foreach(i = 1:length(groups), .combine = list, .multicombine = TRUE) %do%
    {
        group <- groups[i]
        group_comps <- c(group_comps, as.character(group))
        out <- NULL
        for (i in 1:length(pdata[pdata[[group_var]] == group, "sample"])) {
            out <- c(out, paste(pdata[pdata[[group_var]] == group, "sample"][i], pdata[pdata[[group_var]] == "Control", "sample"], sep = "_vs_"),
                     paste(pdata[pdata[[group_var]] == "Control", "sample"], pdata[pdata[[group_var]] == group, "sample"][i], sep = "_vs_"))
        }
        out
    }
    names(res) <- group_comps
    return(res)
}

summarizeBF <- function(group_comps, grouping_var, main_grouping_var, con_bf, pdata, bf_wide_format) {
    if (bf_wide_format) {
        for (i in 1:length(group_comps)) {
            sample_group <- pdata[pdata[[grouping_var]] == names(group_comps)[i], "sample"]
            for (j in 1:length(sample_group)) {
                sample_ID <- sample_group[j]
                tmp <- con_bf[, names(con_bf) %in% group_comps[[i]][grep(pattern = sample_ID, group_comps[[i]])]]
                medians <- apply(tmp, 1, median, na.rm = TRUE)
                assign(as.character(sample_ID),
                       data.frame(con_bf[[main_grouping_var]], medians),
                       envir = .GlobalEnv)
                setnames(get(sample_ID), c(main_grouping_var, sample_ID))
            }
            median_bf <- Reduce(function(...) merge(..., by = main_grouping_var, all = TRUE), lapply(sample_group, get))
            rm(list = sample_group, envir = .GlobalEnv)
            if (i == 1){
                n_sig <- data.table(median_bf[[main_grouping_var]],
                                    n_sig = apply(median_bf[, 2:dim(median_bf)[2], with = F], 1, function(x) length(x[x >= 5 & !is.na(x)])))
                setnames(n_sig, c(main_grouping_var, paste(names(group_comps)[i], "n_sig", sep = "_")))
            } else {
                comps <- data.frame(median_bf[[main_grouping_var]],
                                    n_sig = apply(median_bf[, 2:dim(median_bf)[2], with = F], 1, function(x) length(x[x >= 5 & !is.na(x)])))
                setnames(comps, c(main_grouping_var, paste(names(group_comps)[i], "n_sig", sep = "_")))
                n_sig <- merge(n_sig, comps, by = main_grouping_var, all = TRUE)
            }
        }
    } else {
        if (!"comp" %in% names(con_bf)) {
            con_bf[, comp := paste(sample1, "vs", sample2, sep = "_")]
        }
        for (i in 1:length(group_comps)) {
            sample_group <- pdata[pdata[[grouping_var]] == names(group_comps)[i], "sample"]
            for (j in 1:length(sample_group)) {
                sample_ID <- sample_group[j]
                assign(as.character(sample_ID), con_bf %>%
                           filter(comp %in% group_comps[[i]][grep(pattern = sample_ID, group_comps[[i]])]) %>%
                               group_by_(main_grouping_var) %>%
                                   summarize(n_sig = median(bayes_factor, na.rm = TRUE)),
                       envir = .GlobalEnv)
                setnames(get(sample_ID), c(main_grouping_var, sample_ID))
            }
            median_bf <- Reduce(function(...) merge(..., by = main_grouping_var, all = TRUE), lapply(sample_group, get))
            rm(list = sample_group, envir = .GlobalEnv)
            if (i == 1){
                n_sig <- data.table(median_bf[[main_grouping_var]],
                                    n_sig = apply(median_bf[, 2:dim(median_bf)[2], with = F], 1, function(x) length(x[x >= 5 & !is.na(x)])))
                setnames(n_sig, c(main_grouping_var, paste(names(group_comps)[i], "n_sig", sep = "_")))
            } else {
                comps <- data.table(median_bf[[main_grouping_var]],
                                    n_sig = apply(median_bf[, 2:dim(median_bf)[2], with = F], 1, function(x) length(x[x >= 5 & !is.na(x)])))
                setnames(comps, c(main_grouping_var, paste(names(group_comps)[i], "n_sig", sep = "_")))
                n_sig <- merge(n_sig, comps, by = main_grouping_var, all = TRUE)
            }
        }
    }
    return(n_sig)
}
