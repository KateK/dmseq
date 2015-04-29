##----------------------------------
## required libraries
##----------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(foreach)
library(stringr)

##----------------------------------
## functions
##----------------------------------
find.dm.events <- function(DT) {
    if("gene_symbol" %in% colnames(DT)){
        DT %>%
        filter(abs(delta_psi_mean) >= 0.05,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_n_sig / DM1_n >= 0.25) %>%
                   select(gene_symbol, event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                          DM1_psi_mean, DM1_psi_sd, DM1_n, delta_psi_mean, DM1_n_sig) %>%
                              arrange(desc(abs(delta_psi_mean)))
    } else {
        DT %>%
        filter(abs(delta_psi_mean) >= 0.05,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_n_sig / DM1_n >= 0.25) %>%
                   select(event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                          DM1_psi_mean, DM1_psi_sd, DM1_n, delta_psi_mean, DM1_n_sig) %>%
                              arrange(desc(abs(delta_psi_mean)))
    }
}

##----------------------------------
## load results
##----------------------------------
event_types <- c("nonUTRevents.multi", "ALE", "AFE")

tibialis_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_tibialis_pdata.txt"))
quadricep_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_quadricep_pdata.txt"))
heart_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_heart_pdata.txt"))

nuclear_lamina_genes <- c("LMNA", "")

perturbed_lamina_events <- foreach(i = 1:length(event_types), .combine = rbind, .multicombine = TRUE) %do%
{
    event_type <- event_types[i]
    tibialis_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/tibialis/tibialis", event_type, "results.txt", sep = "_")))
    quadricep_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/quadricep/quadricep", event_type, "results.txt", sep = "_")))
    heart_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/heart/heart", event_type, "results.txt", sep = "_")))
    dm_heart <- find.dm.events(heart_res)
    dm_tibialis <- find.dm.events(tibialis_res)
    dm_quadricep <- find.dm.events(quadricep_res)    
    if (event_type == "nonUTRevents.multi") {
        rbind(dm_heart %>% mutate(Tissue = "Heart", event_type = event_type) %>%
                  filter(gene_symbol == "SYNE1") %>%
                      select(event_name, Tissue, event_type, contains("psi"), contains("_n")),
              dm_quadricep %>% mutate(Tissue = "Quad", event_type = event_type) %>%
                  filter(gene_symbol == "SYNE1") %>%
                      select(event_name, Tissue, event_type, contains("psi"), contains("_n")),
              dm_tibialis %>% mutate(Tissue = "Tibialis", event_type = event_type) %>%
                  filter(gene_symbol == "SYNE1") %>%
                      select(event_name, Tissue, event_type, contains("psi"), contains("_n")))
    } else {
        rbind(dm_heart %>% mutate(Tissue = "Heart", event_type = event_type) %>%
                  filter(event_name == "SYNE1") %>%
                      select(event_name, Tissue, event_type, contains("psi"), contains("_n")),
              dm_quadricep %>% mutate(Tissue = "Quad", event_type = event_type) %>%
                  filter(event_name == "SYNE1") %>%
                      select(event_name, Tissue, event_type, contains("psi"), contains("_n")),
              dm_tibialis %>% mutate(Tissue = "Tibialis", event_type = event_type) %>%
                  filter(event_name == "SYNE1") %>%
                      select(event_name, Tissue, event_type, contains("psi"), contains("_n")))
    }
} 

write.table(perturbed_lamina_events, "~/Projects/DMseq/results/SYNE1_results.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
