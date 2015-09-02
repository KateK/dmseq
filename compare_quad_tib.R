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
        filter(abs(delta_psi_mean) >= 0.05,
               Control_psi_sd < 0.2,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_n_sig / DM1_n >= 0.25) %>%
        select(gene_symbol, event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                      DM1_psi_mean, DM1_psi_sd, DM1_n, delta_psi_mean, DM1_n_sig) %>%
        arrange(desc(abs(delta_psi_mean)))
    } else {
        DT %>%
        filter(abs(delta_psi_mean) >= 0.05,
               Control_psi_sd < 0.2,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_n_sig / DM1_n >= 0.25) %>%
        select(event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                      DM1_psi_mean, DM1_psi_sd, DM1_n, delta_psi_mean, DM1_n_sig) %>%
        arrange(desc(abs(delta_psi_mean)))
    }
}


##----------------------------------
## comparison data
##----------------------------------
event_type <- "nonUTRevents.multi"

allControls_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/allControls/allControls", event_type, "results.txt", sep = "_")))

quad_vs_tibialis <- allControls_res %>%
    select(gene_symbol, event_name, Quad_psi_mean, Quad_n, Tibialis_psi_mean, Tibialis_n, Quad_vs_Tibialis_n_sig, isoforms) %>%
        mutate(delta_psi = abs(Quad_psi_mean - Tibialis_psi_mean)) %>%
            filter(Quad_n / max(Quad_n, na.rm = TRUE) >= 0.75, Tibialis_n / max(Tibialis_n, na.rm = TRUE) >= 0.75,
                   delta_psi > 0.05, Quad_vs_Tibialis_n_sig > 4) %>%
                       arrange(desc(Quad_vs_Tibialis_n_sig))

tibialis_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/tibialis/tibialis", event_type, "results.txt", sep = "_")))
tibialis_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_tibialis_pdata.txt"))

quadricep_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_quadricep_pdata.txt"))
quadricep_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/quadricep/quadricep", event_type, "results.txt", sep = "_")))

heart_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_heart_pdata.txt"))
heart_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/heart/heart", event_type, "results.txt", sep = "_")))

dm_heart <- find.dm.events(heart_res)
dm_tibialis <- find.dm.events(tibialis_res)
dm_quadricep <- find.dm.events(quadricep_res)

dysregulated_events <- Reduce(union, list(dm_tibialis$isoforms, dm_quadricep$isoforms))
gene_set <- intersect(quad_vs_tibialis$isoforms, dysregulated_events)

##----------------------------------
## load event data
##----------------------------------
con_data <- tbl_dt(fread(paste("~/Projects/DMseq/data/allSamples", event_type, "consolidatedSummaries.txt", sep = "_")))
pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_sample_pdata.txt"))
setnames(pdata, c("sample", "diagnosis", "tissue", "group", "read_length"))

con_data <- left_join(con_data, select(pdata, sample, tissue, diagnosis), by = "sample")
con_data <- con_data %>% filter(tissue %in% c("Tibialis", "Quad"))

## filter to keep 'high' confidence psi estimates
f_con_data <- filter(con_data, ci_high - ci_low <= 0.33)

## define sets of all events quatified per tissue and find the intersection of these sets 
tissues <- unique(f_con_data$tissue)

events_by_tissue <- foreach(i = 1:length(tissues), .combine = list, .multicombine = TRUE) %do% {    
    tmp <- f_con_data %>% filter(tissue == tissues[i]) %>% select(isoforms)
    tmp$isoforms
}
names(events_by_tissue) <- tissues
common_events <- Reduce(intersect, events_by_tissue)

## keep only those events quantified in all tissues
common_data <- f_con_data %>% filter(isoforms %in% common_events) %>% select(miso_posterior_mean, tissue, isoforms) %>% as_data_frame

common_data_spread <- tidyr::spread(common_data, key = tissue, value = miso_posterior_mean)

ggplot(common_data_spread) +
    geom_point(aes(x=Tibialis, y=Quad))

