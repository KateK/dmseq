library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

load("~/Projects/DMseq/bin/sum.InformativeMisoCounts.R")

##----------------------------------
## load event data
##----------------------------------
event_type <- "nonUTRevents.multi"

con_data <- tbl_dt(fread(paste("~/Projects/DMseq/data/allSamples", event_type, "consolidatedSummaries.txt", sep = "_")))

alignment_metrics <- tbl_dt(fread("~/Projects/DMseq/data/DM_sample_alignment_metrics.txt"))
pdata <- tbl_dt(fread("~/Projects/DMseq/data/DM_sample_pdata.txt"))
setNames(pdata, c("sample", "diagnosis", "tissue", "group", "read_length"))

con_data <- con_data %>% mutate(informativeCounts = sum.InformativeMisoCounts(counts))

head(con_data,1000) %>%
    group_by(sample) %>%
    filter(informativeCounts >= 20) %>%
    summarize(n_events = n())

n_events <- left_join(con_data %>% group_by(sample) %>% filter(informativeCounts >= 20) %>% summarize(n_events_covered = n()),
                      select(alignment_metrics, sample, fastq_read_counts, concordant_unique_read_counts, annotated_junction_read_counts),
                      by = "sample")

n_events <- left_join(n_events, pdata, by = "sample")

ggplot(n_events, aes(y = n_events_covered, x = fastq_read_counts)) +
    geom_point(aes(colour = factor(tissue))) +
        geom_smooth(colour = "black") +
            labs(x = "# of reads from sequencing",
                 y = "# of splicing events covered")
