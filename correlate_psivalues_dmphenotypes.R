##----------------------------------
## sample metadata
##----------------------------------
pdata_file <- "~/Analysis_Projects/mbnl1_dose_paper/data/DM_tibialis_pdata.txt"

##----------------------------------
## nonUTRevents.multi
##----------------------------------
consolidated_summary_file <- "~/MISO/single_end_mode/summaries/tibialis/nonUTRevents.multi_consolidatedSummaries.txt"
event_key_file <- "~/alt.splicing.annotations/hg19_Aug2014/nonUTRevents.multi.hg19.to.ensGene.txt"

cor_data <- calc_correlations(consolidated_summary_file, pdata_file,
                              gff_file = NULL, event_key_file = event_key_file,
                              grouping_var = "isoforms", counts_filter = TRUE)

write.table(cor_data, file = "~/Analysis_Projects/DMseq/results/nonUTRevents.multi_cor.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
## ALE
##----------------------------------
consolidated_summary_file <- "~/MISO/single_end_mode/summaries/tibialis/ALE_consolidatedSummaries.txt"
event_key_file <- NULL

cor_data <- calc_correlations(consolidated_summary_file, pdata_file,
                              gff_file = NULL, event_key_file = event_key_file,
                              grouping_var = "isoforms", counts_filter = TRUE)

write.table(cor_data, file = "~/Analysis_Projects/DMseq/results/ALE_cor.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
## CLCN1_multi
##----------------------------------
consolidated_summary_file <- "~/MISO/single_end_mode/summaries/tibialis/CLCN1_multi_consolidatedSummaries.txt"
event_key_file <- NULL


cor_data <- calc_correlations(consolidated_summary_file, pdata_file,
                              gff_file = NULL, event_key_file = event_key_file,
                              grouping_var = "isoforms", counts_filter = TRUE)

write.table(cor_data, file = "~/Analysis_Projects/DMseq/results/CLCN1_multi_cor.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##----------------------------------
##----------------------------------
## Analysis
##----------------------------------
##----------------------------------
calc_correlations <- function(consolidated_summary_file, pdata_file,
                              gff_file = NULL, event_key_file = NULL,
                              grouping_var = "event_name", counts_filter = TRUE)
{
    require(data.table)
    require(dplyr)
    
    con_summaries <- tbl_dt(fread(consolidated_summary_file))
    con_summaries[, ci_width := ci_high - ci_low]
    con_summaries[, informativeCounts := sum.InformativeMisoCounts(con_summaries$counts)]
    
    pdata <- read.table(pdata_file, header = TRUE, row.names = NULL, sep = "\t") 
    names(pdata)[1] <- "sample"

    if(counts_filter) {
        con_summaries <- con_summaries[informativeCounts >= 20,]        
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
        named.events <- con_summaries[ ,c("event_name", "isoforms"), with = F] %>% unique()
    }
    
    ## combine consolidated data with pdata ##
    consolidated_data <- tbl_dt(merge(con_summaries, pdata, by = "sample"))

    ## Correlate psi values with strength measurements ##
    con_data <- consolidated_data %>% group_by_(grouping_var) %>% mutate(n = sum(!is.na(miso_posterior_mean)))
    cor_data_strength <- con_data %>%
        filter(n > 10) %>%
            group_by_(grouping_var) %>%
                summarize(cor_HG.QMT = cor.test(miso_posterior_mean, HG.QMT, method = "pearson", use = "pairwise.complete.obs")$estimate,
                          pval_HG.QMT = cor.test(miso_posterior_mean, HG.QMT, method = "pearson", use = "pairwise.complete.obs")$p.value,
                          n_complete_HG.QMT = sum(complete.cases(miso_posterior_mean, HG.QMT)),
                          cor_ADF.QMT = cor.test(miso_posterior_mean, ADF.QMT, method = "pearson", use = "pairwise.complete.obs")$estimate,
                          pval_ADF.QMT = cor.test(miso_posterior_mean, ADF.QMT, method = "pearson", use = "pairwise.complete.obs")$p.value,
                          n_complete_ADF.QMT = sum(complete.cases(miso_posterior_mean, ADF.QMT)),
                          cor_Actual.Strength.6pt.Scale = cor.test(miso_posterior_mean, Actual.Strength.6pt.Scale, method = "spearman", use = "pairwise.complete.obs")$estimate,
                          pval_Actual.Strength.6pt.Scale = cor.test(miso_posterior_mean, Actual.Strength.6pt.Scale, method = "spearman", use = "pairwise.complete.obs")$p.value,
                          n_complete_Actual.Strength.6pt.Scale = sum(complete.cases(miso_posterior_mean, Actual.Strength.6pt.Scale)))

    ## Correlate psi values with repeat length ##
    con_data <- consolidated_data %>% filter(!is.na(as.numeric(as.character(n.CTG.repeats)))) %>% group_by_(grouping_var) %>% mutate(n = sum(!is.na(miso_posterior_mean)))
    cor_data_repeats <- con_data %>%
        group_by_(grouping_var) %>%
            filter(n > 4) %>%
                summarize(cor_n.CTG.repeats = cor.test(miso_posterior_mean, as.numeric(as.character(n.CTG.repeats)), method = "pearson", use = "pairwise.complete.obs")$estimate,
                          pval_n.CTG.repeats = cor.test(miso_posterior_mean, as.numeric(as.character(n.CTG.repeats)), method = "pearson", use = "pairwise.complete.obs")$p.value,
                          n_complete_n.CTG.repeats = sum(complete.cases(miso_posterior_mean, as.numeric(as.character(n.CTG.repeats)))))
    ## Merge results ##
    cor_data <- merge(cor_data_strength, cor_data_repeats, by = grouping_var, all = TRUE)
    if(exists("named.events")) {
        cor_data <- merge(cor_data, named.events, by = grouping_var)
    }
    if(exists("event.key")) {
        cor_data <- merge(cor_data, event.key, by = "event_name", all.x = TRUE)
    }
    return(cor_data)
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
