##----------------------------------
## required libraries
##----------------------------------
library(dplyr)
library(data.table)
library(reshape2)
library(foreach)

pdata_file <- "~/Projects/DM_tibialis_analysis/data/DM_tibialis_pdata.txt"

##----------------------------------
## nonUTRevents.multi
##----------------------------------
consolidated_summary_file <- "~/Projects/DM_tibialis_analysis/data/nonUTRevents.multi_consolidatedSummaries.txt"
event_key_file <- "~/Projects/DM_tibialis_analysis/data/nonUTRevents.multi.hg19.to.ensGene.txt"

##----------------------------------
## ALE
##----------------------------------
consolidated_summary_file <- "~/Projects/DM_tibialis_analysis/data/ALE_consolidatedSummaries.txt"
event_key_file <- NULL

##----------------------------------
##----------------------------------
## Analysis
##----------------------------------
##----------------------------------

### Load files into memory ###
con_summaries <- tbl_dt(fread(consolidated_summary_file))
con_summaries[, ci_width := ci_high - ci_low]

pdata <- read.table(pdata_file, header = TRUE, row.names = NULL, sep = "\t") 
names(pdata)[1] <- "sample"

### Map MISO annotations to gene symbols ###
if (!is.null(gff_file)) {
    gff <- fread(gff_file)
    gff <- gff %>% dplyr::filter(V3 == "gene")
    split <- strsplit(as.character(gff$V9), "[;=]")
    named.events <- data.table(do.call(rbind.data.frame, lapply(split, function(x) { 
        x[c(2,8,10,12)]
    })))
}

if (!is.null(event_key_file)) {
    named.events <- fread(event_key_file)
    setnames(named.events, c("event_name", "ensg_id", "gene_symbol", "desc"))
}

### combine consolidated data with pdata ###
consolidated_data <- tbl_dt(merge(con_summaries, pdata, by="sample"))
consolidated_data$n.CTG.repeats <- as.numeric(as.character(consolidated_data$n.CTG.repeats))

### aov w/ gender ###
grouping_var <- "isoforms"

con_data_wide <- dcast(consolidated_data, sample + Gender + Age + n.CTG.repeats ~ isoforms, value.var = "miso_posterior_mean")

### only consider events with most (50) of the samples having data
filtered_events <- which(colSums(!is.na(con_data_wide[,1:dim(con_data_wide)[2]]), na.rm = T) >= 50)
con_data_wide <- con_data_wide[, filtered_events]

##----------------------------------
## analysis of gender affects on psi
##----------------------------------
aov_res <- foreach(i = 4:dim(con_data_wide)[2], .combine = c, .multicombine = TRUE) %do%
{
    event <- names(con_data_wide)[i]
    anova(aov(get(event) ~ Gender, data = con_data_wide))[["Pr(>F)"]][1]
}

aov_res <- setNames(aov_res, names(con_data_wide)[4:length(names(con_data_wide))])

## Benjamini & Hochberg correction
aov_res_adj <- p.adjust(aov_res, method = "BH")
## Bonferroni correction
aov_res_adj <- p.adjust(aov_res, method = "bonferroni")

##----------------------------------
## analysis of age affects on psi
##----------------------------------
lm_res <- foreach(i = 4:dim(con_data_wide)[2], .combine = c, .multicombine = TRUE) %do%
{
    event <- names(con_data_wide)[i]
    summary(lm(get(event) ~ Age, data = con_data_wide))[["coefficients"]][2,4]
}

lm_res <- setNames(lm_res, names(con_data_wide)[4:length(names(con_data_wide))])

## Benjamini & Hochberg correction
lm_res_adj <- p.adjust(lm_res, method = "BH")
## Bonferroni correction
lm_res_adj <- p.adjust(aov_res, method = "bonferroni")
