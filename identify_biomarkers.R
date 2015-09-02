##----------------------------------
## required libraries
##----------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(foreach)
library(ggplot2)
library(RColorBrewer)

##----------------------------------
## functions
##----------------------------------
find.dm.events <- function(DT) {
    if("gene_symbol" %in% colnames(DT)){
        DT %>%
        filter(abs(Control_psi_mean - DM1_psi_mean) >= 0.05,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_vs_Control_n_sig / DM1_n >= 0.25) %>%
        select(gene_symbol, event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                      DM1_psi_mean, DM1_psi_sd, DM1_n, DM1_vs_Control_n_sig) %>%
        arrange(desc(abs(Control_psi_mean - DM1_psi_mean)))
    } else {
        DT %>%
        filter(abs(Control_psi_mean - DM1_psi_mean) >= 0.05,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_vs_Control_n_sig / DM1_n >= 0.25) %>%
        select(event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                      DM1_psi_mean, DM1_psi_sd, DM1_n, DM1_vs_Control_n_sig) %>%
        arrange(desc(abs(Control_psi_mean - DM1_psi_mean)))
    }
}


##----------------------------------
## input files
##----------------------------------
event_type <- "nonUTRevents.multi"
tissue <- "tibialis"

pdata_f <- paste("~/Projects/DMseq/data/DM_", tissue, "_full_pdata.txt", sep = "")
con_f <- paste("~/Projects/DMseq/data/", tissue, "/", tissue, "_", event_type, "_consolidatedSummaries.txt", sep = "")
res_f <- paste("~/Projects/DMseq/results/", tissue, "/", tissue, "_", event_type, "_results.txt", sep = "")
cor_f <- paste("~/Projects/DMseq/results/", tissue, "/", tissue, "_", event_type, "_cor.txt", sep = "")

##----------------------------------
## Analysis
##----------------------------------
con_data <- fread(con_f) %>% tbl_dt()

res_data <- fread(res_f) %>% tbl_dt()
dm_events <- find.dm.events(res_data)

cor_data <- fread(cor_f) %>% tbl_dt()
pdata <- fread(pdata_f) %>% tbl_dt()

cor.test(pdata$HG.QMT, pdata$ADF.QMT, use = "na.or.complete")

ggplot(pdata, aes(x = HG.QMT, ADF.QMT)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "longdash")


cor_data %>%
    filter(isoforms %in% dm_events$isoforms) %>%
    select(event_name, gene_symbol, contains("cor")) %>%
    arrange(desc(abs(cor_ADF.QMT))) %>% head(20)
