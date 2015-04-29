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

sample_pdata <- tbl_dt(fread("~/Analysis_Projects/DMseq/data/DM_sample_pdata.txt"))
allControls_res <- tbl_dt(fread(paste("~/Analysis_Projects/DMseq/results/allControls/allControls", event_type, "results.txt", sep = "_")))

quad_vs_tibialis <- allControls_res %>%
    select(gene_symbol, event_name, Quad_psi_mean, Quad_n, Tibialis_psi_mean, Tibialis_n, Quad_vs_Tibialis_n_sig, isoforms) %>%
        mutate(delta_psi = abs(Quad_psi_mean - Tibialis_psi_mean)) %>%
            filter(Quad_n / max(Quad_n, na.rm = TRUE) >= 0.75, Tibialis_n / max(Tibialis_n, na.rm = TRUE) >= 0.75,
                   delta_psi > 0.05, Quad_vs_Tibialis_n_sig > 4) %>%
                       arrange(desc(Quad_vs_Tibialis_n_sig))

tibialis_res <- tbl_dt(fread(paste("~/Analysis_Projects/DMseq/results/tibialis/tibialis", event_type, "results.txt", sep = "_")))
tibialis_pdata <- tbl_dt(fread("~/Analysis_Projects/DMseq/data/tibialis/DM_tibialis_pdata.txt"))

quadricep_pdata <- tbl_dt(fread("~/Analysis_Projects/DMseq/data/quadricep/DM_quadricep_pdata.txt"))
quadricep_res <- tbl_dt(fread(paste("~/Analysis_Projects/DMseq/results/quadricep/quadricep", event_type, "results.txt", sep = "_")))

heart_pdata <- tbl_dt(fread("~/Analysis_Projects/DMseq/data/heart/DM_heart_pdata.txt"))
heart_res <- tbl_dt(fread(paste("~/Analysis_Projects/DMseq/results/heart/heart", event_type, "results.txt", sep = "_")))

dm_heart <- find.dm.events(heart_res)
dm_tibialis <- find.dm.events(tibialis_res)
dm_quadricep <- find.dm.events(quadricep_res)

dysregulated_events <- Reduce(union, list(dm_tibialis$isoforms, dm_quadricep$isoforms))
gene_set <- intersect(quad_vs_tibialis$isoforms, dysregulated_events)

##----------------------------------
## nonUTRevents.nonmulti
##----------------------------------
con_f <- "~/MISO/single_end_mode/summaries/tibialis/nonUTRevents.nonmulti_consolidatedSummaries.txt"
res_f <- "~/Analysis_Projects/DMseq/results/nonUTRevents.nonmulti_results.txt"
cor_f <- "~/Analysis_Projects/DMseq/results/nonUTRevents.nonmulti_strength_cor.txt"
repeat_length_cor_f <- "~/Analysis_Projects/DMseq/results/nonUTRevents.nonmulti_repeat_length_cor.txt"
event_key_file <- "~/alt.splicing.annotations/hg19_Aug2014/nonUTRevents.nonmulti.hg19.to.ensGene.txt"

##----------------------------------
## nonUTRevents.multi
##----------------------------------
con_f <- "~/MISO/single_end_mode/summaries/tibialis/nonUTRevents.multi_consolidatedSummaries.txt"
res_f <- "~/Analysis_Projects/DMseq/results/nonUTRevents.multi_results.txt"
cor_f <- "~/Analysis_Projects/DMseq/results/nonUTRevents.multi_cor.txt"
event_key_file <- "~/alt.splicing.annotations/hg19_Aug2014/nonUTRevents.multi.hg19.to.ensGene.txt"


##----------------------------------
## ALE
##----------------------------------
con_f <- "~/MISO/single_end_mode/summaries/tibialis/ALE_consolidatedSummaries.txt"
res_f <- "~/Analysis_Projects/DMseq/results/ALE_results.txt"
cor_f <- "~/Analysis_Projects/DMseq/results/ALE_strength_cor.txt"

res_data <- fread(res_f)
cor_data <- fread(cor_f)
con_data <- fread(con_f)

events <- fread(event_key_file)
setnames(events, c("event_name","ensg_id", "gene_symbol", "desc"))

pdata <- read.table(pdata_file, header = TRUE, row.names = NULL, sep = "\t") 
names(pdata)[1] <- "sample"

con_data <- merge(con_data, pdata, by = "sample")

tmp <- cor_data %>% filter(n_complete_HG.QMT > 25, n_complete_ADF.QMT > 25, n_complete_Actual.Strength.6pt.Scale > 25, pval_HG.QMT < 0.05, pval_ADF.QMT < 0.05, pval_Actual.Strength.6pt.Scale < 0.05)
tmp <- cor_data %>% filter(n_complete_HG.QMT > 25, n_complete_ADF.QMT > 25, n_complete_Actual.Strength.6pt.Scale > 25)
tmp <- tmp %>% arrange(abs(cor_Actual.Strength.6pt.Scale), abs(cor_HG.QMT), abs(cor_ADF.QMT))
tmp <- merge(tmp, events, by = "event_name")

arrayevents <- read.table("~/Projects/DM_tibialis_analysis/data/biomarker_event_key.txt", header = TRUE, sep = "\t")

event <- "chr19:7152737:7152938:-@chr19:7150508:7150543:-@chr19:7142827:7143101:-"
ggplot(con_data[event_name == event,], aes(x = Actual.Strength.6pt.Scale.Norm7, y = miso_posterior_mean)) + geom_point() + ylab("PSI") + xlab("MRC scale") + ggtitle("TTN") + theme_grey(30)
ggplot(con_data[event_name == event,], aes(x = as.numeric(as.character(n.CTG.repeats)), y = miso_posterior_mean)) + geom_point() + ylab("PSI") + xlab("# of repeats") + theme_grey(30)

### must set up what feild is being shuffled manually right now
shuffle_and_calc_cors <- function(con_data, nshuffles) {
    x <- con_data %>% group_by(event_name) %>% mutate(n = sum(!is.na(miso_posterior_mean)))
    tmp <- foreach(i = 1:nshuffles, .combine = cbind) %do%
    {
        x$shuffled <- sample(x$Actual.Strength.6pt.Scale, length(x$Actual.Strength.6pt.Scale), replace = FALSE)
        x %>% filter(n > 3) %>% group_by(event_name) %>% summarize(cor = cor(miso_posterior_mean, shuffled, method = "pearson", use = "pairwise.complete.obs")) %>% select(cor)
    }
    apply(tmp, 1, function(x) mean(x, na.rm = TRUE))
}

shuffled_MRC_cors <- shuffle_and_calc_cors(con_data, 100)
sampled_shuf_cors <- sample(shuffled_MRC_cors, length(cor_data$cor_Actual.Strength.6pt.Scale), replace = FALSE)

shuffle_and_calc_cors <- function(con_data, nshuffles) {
    x <- con_data %>% group_by(event_name) %>% mutate(n = sum(!is.na(miso_posterior_mean)))
    x <- x %>% filter(!is.na(n.CTG.repeats))
    tmp <- foreach(i = 1:nshuffles, .combine = cbind) %do%
    {
        x$shuffled <- sample(as.numeric(as.character(x$n.CTG.repeats)), length(x$n.CTG.repeats), replace = FALSE)
        x %>% filter(n > 3) %>% group_by(event_name) %>% summarize(cor = cor(miso_posterior_mean, shuffled, method = "pearson", use = "pairwise.complete.obs")) %>% select(cor)
    }
    apply(tmp, 1, function(x) mean(x, na.rm = TRUE))
}


shuffled_repeat_cors <-  shuffle_and_calc_cors(con_data, 100)
sampled_shuf_cors <- sample(shuffled_repeat_cors, length(rep_data$cor), replace = FALSE)

### plot observed cor values versus values dervided from shuffled data
ggplot() + geom_freqpoly(binwidth = 0.01, colour = 'black', aes(x = abs(sampled_shuf_cors))) +
    geom_freqpoly(binwidth = 0.01, colour = 'red', aes(x = abs(cor_data$cor_Actual.Strength.6pt.Scale))) +
    xlim(0,1) + xlab("Pearson's r") + theme_grey(30)

ggplot() + geom_freqpoly(binwidth = 0.01, colour = 'black', aes(x = abs(shuffled_cors))) + xlim(0,1)
ggplot() + geom_freqpoly(binwidth = 0.01, colour = 'red', aes(x = abs(cor_data$cor_Actual.Strength.6pt.Scale))) + xlim(0,1)

###
psum <- function(..., na.rm=FALSE) {
    x <- list(...)
    rowSums(matrix(unlist(x), ncol=length(x)), na.rm=na.rm)
}

### ci_width doesn't correlate well with number of informative reads
con_data[, ci_width := ci_high - ci_low]
con_data[, inform_reads := psum(reads_first_iso, reads_second_iso, reads_both_iso, na.rm = TRUE)]

cor.test(con_data$ci_width, con_data$inform_reads, use = "pairwise.complete.obs")

### count sig events for each event_type
tmp <- res_data %>% select(event_name, gene_symbol, Control_psi_mean, Control_ci_width_mean, Control_n, DM1_psi_mean, DM1_ci_width_mean, DM1_n, delta_psi_mean, DM1_n_sig) %>% filter(DM1_n_sig > 0)
non_SS <- tmp[grep("[|]", tmp$event_name, invert = TRUE)]

alt_SS <- tmp[grep("[|]", tmp$event_name)]
SE <- tmp[grep("^chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]$", tmp$event_name)]
MXE <- tmp[grep("^chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]$", tmp$event_name)]
RI <- tmp[grep("^chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][-+]@chr[0-9XYM]+[:][0-9]+[:][0-9]+[:][+-]$", tmp$event_name)]

dim(filter(alt_SS, abs(delta_psi_mean) >= 0.05, DM1_n > 23, Control_n > 6, DM1_n_sig / (DM1_n * Control_n) >=  0.33))
dim(filter(SE, abs(delta_psi_mean) >= 0.05, DM1_n > 23, Control_n > 6, DM1_n_sig / (DM1_n * Control_n) >=  0.33))
dim(filter(MXE, abs(delta_psi_mean) >= 0.05, DM1_n > 23, Control_n > 6, DM1_n_sig / (DM1_n * Control_n) >=  0.33))
dim(filter(RI, abs(delta_psi_mean) >= 0.05, DM1_n > 23, Control_n > 6, DM1_n_sig / (DM1_n * Control_n) >=  0.33))

###
res_sig <- filter(res_data, abs(delta_psi_mean) >= 0.05, DM1_n > 23, Control_n > 6, DM1_n_sig / (DM1_n * Control_n) >=  0.33)

shuffle_and_calc_cors <- function(con_data, shuffle_var, nshuffles) {
    n <- con_data %>% filter(event_name %in% res_sig$event_name) %>% group_by(event_name) %>% tally()
    x <-  merge(con_data, n, by = "event_name")
    tmp <- foreach(i = 1:nshuffles, .combine = cbind) %do%
    {
        x$shuffled <- sample(x[[shuffle_var]], length(x[[shuffle_var]]), replace = FALSE)
        x %>% filter(n > 3) %>% group_by(event_name) %>% summarize(cor = cor(miso_posterior_mean, shuffled, method = "pearson", use = "pairwise.complete.obs")) %>% select(cor)
    }
    apply(tmp, 1, function(x) mean(x, na.rm = TRUE))
}

shuffled_MRC_cors <- shuffle_and_calc_cors(con_data, "Actual.Strength.6pt.Scale", 10)
actual_cors <- cor_data %>% filter(event_name %in% res_sig$event_name)
actual_cors <- actual_cors$cor_Actual.Strength.6pt.Scale
sampled_shuf_cors <- sample(shuffled_MRC_cors, length(actual_cors), replace = FALSE)

ggplot() + geom_freqpoly(binwidth = 0.1, colour = 'black', aes(x = abs(sampled_shuf_cors))) +
    geom_freqpoly(binwidth = 0.1, colour = 'red', aes(x = abs(actual_cors))) +
    xlim(-0.05,1) + xlab("Pearson's r") + theme_grey(30)

ggplot() + geom_freqpoly(binwidth = 0.1, colour = 'black', aes(x = sampled_shuf_cors)) +
    geom_freqpoly(binwidth = 0.1, colour = 'red', aes(x = actual_cors)) +
    xlab("Pearson's r") + theme_grey(30)

ggplot() + geom_histogram(binwidth = 0.01, fill = 'black', alpha = 0.2, aes(x = abs(sampled_shuf_cors))) +
    geom_histogram(binwidth = 0.01, fill = 'red', alpha = 0.2, aes(x = abs(actual_cors))) +
    xlab("Pearson's r") + theme_grey(30)

ggplot() + geom_histogram(binwidth = 0.01, fill = 'black', alpha = 0.2, aes(x = sampled_shuf_cors)) +
    geom_histogram(binwidth = 0.01, fill = 'red', alpha = 0.2, aes(x = actual_cors)) +
    xlab("Pearson's r") + theme_grey(30)

