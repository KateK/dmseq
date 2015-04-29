##----------------------------------
## required libraries
##----------------------------------
library(foreach)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

library(RColorBrewer)
library(gplots)

##----------------------------------
## functions
##----------------------------------
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


tissue.color.map <- function(tissues)
{
    color_key$color[match(tissues, color_key$tissue)]
}

diagnosis.color.map <- function(diagnosis, colors)
{
    if(diagnosis == "DM1") colors[1]
    else if(diagnosis == "Control") colors[2]
    else if(diagnosis == "Proto_DM1") colors[3]
}

find.dm.events <- function(DT) {
    DT %>%
        filter(abs(delta_psi_mean) >= 0.05,
               Control_psi_sd < 0.2,
               Control_n / max(Control_n, na.rm = TRUE)  >= 0.75,
               DM1_n / max(DM1_n, na.rm = TRUE) >= 0.75,
               DM1_n_sig / DM1_n >= 0.25) %>%
        select(gene_symbol, event_name, isoforms, Control_psi_mean, Control_psi_sd, Control_n,
                      DM1_psi_mean, DM1_psi_sd, DM1_n, delta_psi_mean, DM1_n_sig) %>%
        arrange(desc(abs(delta_psi_mean)))
}

##----------------------------------
## load event data
##----------------------------------
event_type <- "nonUTRevents.multi"
event_type <- "ALE"
event_type <- "AFE"

con_data <- tbl_dt(fread(paste("~/Projects/DMseq/data/allSamples", event_type, "consolidatedSummaries.txt", sep = "_")))
pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_sample_pdata.txt"))
setNames(pdata, c("sample", "diagnosis", "tissue", "group", "read_length"))

con_data <- left_join(con_data, select(pdata, sample, tissue, diagnosis), by = "sample")
con_data <- con_data %>% mutate(ci_width = ci_high - ci_low)

## filter to keep 'high' confidence psi estimates
f_con_data <- filter(con_data, ci_width < 0.33)

## define sets of all events quatified per tissue and find the intersection of these sets 
tissues <- unique(pdata$tissue)
events_by_tissue <- foreach(i = 1:length(tissues), .combine = list, .multicombine = TRUE) %do% {    
    tmp <- f_con_data %>% filter(tissue == tissues[i]) %>% select(isoforms)
    tmp$isoforms
}
names(events_by_tissue) <- tissues
common_events <- Reduce(intersect, events_by_tissue)

## keep only those events quantified in all tissues
common_data <- f_con_data %>% filter(isoforms %in% common_events)
## common_data[, informativeCounts := sum.InformativeMisoCounts(common_data$counts)]

##--------------------------------
## Heatmap colors
##--------------------------------
my_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(1000)
## my_palette <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(1000)

## define colors to identify tissues
color_key <- data.table(tissue = tissues, color = rainbow(length(tissues)))

##----------------------------------
## results
##----------------------------------
tibialis_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/tibialis/tibialis", event_type, "results.txt", sep = "_")))
tibialis_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_tibialis_pdata.txt"))

quadricep_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_quadricep_pdata.txt"))
quadricep_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/quadricep/quadricep", event_type, "results.txt", sep = "_")))

heart_pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_heart_pdata.txt"))
heart_res <- tbl_dt(fread(paste("~/Projects/DMseq/results/heart/heart", event_type, "results.txt", sep = "_")))

dm_heart <- find.dm.events(heart_res)
dm_tibialis <- find.dm.events(tibialis_res)
dm_quadricep <- find.dm.events(quadricep_res)

##--------------------------------
## PSI heatmap - nonUTRevents
##--------------------------------
## randomly select splicing events
## gene_set <- sample(unique(common_data$isoforms), 2000, replace = FALSE)

## select events that show evidence of dysregulation in the quad, tibialis or heart
dysregulated_events <- Reduce(union, list(dm_heart$isoforms, dm_tibialis$isoforms, dm_quadricep$isoforms))
gene_set <- intersect(common_events, dysregulated_events)

f_common_data <- common_data %>% filter(isoforms %in% gene_set) %>% select(sample, isoforms, miso_posterior_mean, counts)
## filter out points w/o sufficient coverage
f_common_data[, informativeCounts := sum.InformativeMisoCounts(f_common_data$counts)]
f_common_data <- f_common_data %>% filter(informativeCounts >= 20) %>% select(sample, isoforms, miso_posterior_mean)

## create event X sample matrix containing psi values
psi_table <- spread(f_common_data, key = sample, value = miso_posterior_mean)

## plottting
sample_tissues <- pdata[match(names(select(psi_table, -isoforms)), sample)]$tissue
sample_diagnosis <- pdata[match(names(select(psi_table, -isoforms)), sample)]$diagnosis
tissue_colors <- tissue.color.map(sample_tissues)

postscript(file = paste("~/Projects/DMseq/results/", event_type, "_heatmap.eps", sep = ""))

heatmap.2(as.matrix(select(psi_table, -isoforms)),
          labCol = paste(names(select(psi_table, -isoforms)), sample_diagnosis, sample_tissues, sep = "-"),
          labRow = "",
          key = TRUE, trace = "none", density.info = "none", keysize = 1,
          key.xlab = "PSI", key.title = NA,
          dendrogram = "both", Colv = TRUE, Rowv = TRUE,
          col = my_palette,
          ColSideColors = tissue_colors,
          margins = c(12, 8))

dev.off()

##--------------------------------
## Cor heatmap - all nonUTRevents
##--------------------------------
## cor for all events
## psi_table_bg <- spread(select(common_data, sample, isoforms, miso_posterior_mean), key = sample, value = miso_posterior_mean)
## sample_cors <- cor(select(psi_table_bg, -isoforms), use = "pairwise.complete.obs", method = "spearman")

## cor of events that show evidence of dysregulation in the quad, tibialis or heart
sample_cors <- cor(select(psi_table, -isoforms), use = "pairwise.complete.obs", method = "spearman")

## plotting 
postscript(file = paste("~/Projects/DMseq/results/sample_cor_", event_type, "_heatmap.eps", sep = ""))

heatmap.2(as.matrix(sample_cors),
          labCol = paste(names(select(psi_table_bg, -isoforms)), sample_diagnosis, sample_tissues, sep = "-"),
          labRow = "",
          key = TRUE, trace = "none", density.info = "none", keysize = 1,
          key.xlab = "Cor", key.title = NA,
          dendrogram = "both", Colv = TRUE, Rowv = TRUE,
          col = my_palette,
          ColSideColors = tissue_colors,
          RowSideColors = tissue_colors,
          margins = c(12, 8))

dev.off()
