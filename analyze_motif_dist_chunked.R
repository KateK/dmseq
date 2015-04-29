##----------------------------------
## required libraries
##----------------------------------
library(foreach)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggbio)
## library(TxDb.Hsapiens.UCSC.hg19.knownGene)

##----------------------------------
## sequence data
##----------------------------------
sequences <- tbl_dt(fread("~/Projects/mbnl1_dose_paper/data/heatmap_gene_seqs_400bp.txt"))
sequences <- tbl_dt(gather(sequences, 2:4, key = "Region", value = "Seq"))

bg_sequences <- tbl_dt(fread("~/Projects/mbnl1_dose_paper/data/background_gene_seqs_400bp.txt"))

bg_sequences <- tbl_dt(tidyr::gather(bg_sequences, 2:4, key = "Region", value = "Seq"))
bg_sequences <- filter(bg_sequences, !Gene_Symbol %in% c(sequences$Gene_Symbol, NA))
all_sequences <- rbind(sequences, bg_sequences)

events_bed <- tbl_dt(fread("~/Projects/mbnl1_dose_paper/data/heatmap_gene_bed_400bp.txt"))

bg_events_bed <- tbl_dt(fread("~/Projects/mbnl1_dose_paper/data/background_gene_bed_400bp.txt"))
bg_events_bed <- filter(bg_events_bed, !name %in% c(events_bed$name, NA))

events_summary <- tbl_dt(fread("~/Projects/mbnl1_dose_paper/results/heatmap_genes_summary.txt"))
events_summary$reg_type <- unlist(lapply(events_summary$delta_psi_mean, function(x) if (sign(x) == 1) "Negative" else "Positive"))

##----------------------------------
## define event list by type
##----------------------------------
events <- unique(all_events_bed$name)
pos_events <-events_summary[reg_type == "Positive",]$gene_symbol
neg_events <-events_summary[reg_type == "Negative",]$gene_symbol
bg_events <- unique(bg_events_bed$name)

##----------------------------------
## mostif distribution analysis 
##----------------------------------
## Bind-n-Seq 4-mers
mbnl_motifs <- paste(c("GCTT", "CGCT", "TGCT", "GCGC", "CCGC", "CTGC", "GCTA", "ACGC", "CGCA", "AGCT", "TTGC", "CAGC"), collapse = "|")
celf_motifs <- paste(c("TGTT", "ATGT", "TTGT", "TGTC", "GTGT", "TGTA", "GTTT", "TGTG", "GTCT", "TTTT"), collapse = "|")

## YGCY 
mbnl_motifs <- paste(c("TGCC", "CGCT", "TGCT", "CGCC"), collapse = "|")

min_motif_size <- 4

## sliding window over sequence; counts motifs per chunk
intron_chunk_size <- 10

motifs_per_event <- foreach(i = 1:length(events), .combine = list, .multicombine = TRUE, .maxcombine = length(events) + 1) %do%
{
    event <- all_sequences %>% filter(Gene_Symbol == events[i])
    regions <- event$Region
    motifs_per_region <- foreach(j = 1:length(regions), .combine = list, .multicombine = TRUE) %do%
    {
        region <- event %>% filter(Region == regions[j])
        if (regions[j] == "EX_Seq") {
            chunk_size <- max(round(nchar(region$Seq) / 4), min_motif_size)
            starts <- seq(1, nchar(region$Seq) - chunk_size, by = chunk_size)
        } else {
            chunk_size <- intron_chunk_size
            starts <- seq(1, nchar(region$Seq) - chunk_size, by = 1)
        }
        ## starts <- seq(1, nchar(region$Seq) - chunk_size, by = 1)
        n <- length(starts)
        motifs_per_chunk <- foreach(k = 1:n, .combine = c, .multicombine = TRUE, .maxcombine = n + 1, .verbose = FALSE) %do%
        { 
            start <- starts[k]
            end <- min(start + chunk_size - 1, nchar(region$Seq))
            chunk <- substr(region$Seq, start, end)
            sum(str_count(mbnl_motifs, chunk))
            setNames(sum(str_count(chunk, mbnl_motifs)), start)
        }
    }
    motifs_per_region <- setNames(motifs_per_region, regions)
    motifs_per_region
}

motifs_per_event <- setNames(motifs_per_event, events)

##----------------------------------
## transform data
##----------------------------------
all_events <- foreach(i=1:length(events), .combine = cbind, .multicombine = FALSE) %do%
{
    event <- events[i]
    event_info <- events_bed %>% filter(name == event)    
    if (i == 1) {
        ## shift upstream positions -20 & set to be negative
        ## shift downstream positions +20
        ## model exon from pos -20 to +20
        event_dat <- data_frame(start = c(-as.numeric(names(motifs_per_event[[event]]$UI_Seq)) - 20, as.numeric(names(motifs_per_event[[event]]$DI_Seq)) + 20),
                                n_motifs = as.numeric(c(motifs_per_event[[event]]$UI_Seq, motifs_per_event[[event]]$DI_Seq)))
        names(event_dat) <- c("start", event)
    } else {
        event_dat <- data_frame(n_motifs = as.numeric(c(motifs_per_event[[event]]$UI_Seq, motifs_per_event[[event]]$DI_Seq)))
        names(event_dat) <- c(event)
    }
    event_dat
}

all_events <- tbl_df(all_events)
all_events <- all_events %>% mutate(end = start + 1, chr = "chr0")

##----------------------------------
## plot with ggbio
##----------------------------------
## n motifs NOT averaged across position
gr_model <- GRanges("chr0", IRanges(start = c(-500, -20, 500), width = c(50, 40, 50)))
p_model <- autoplot(gr_model, geom = "alignment", main.geom = "arrowrect", gap.geom = "segment")

selection <- dim(all_events)[2] - 2
all_events_long <- gather(all_events, 2:selection , key = "gene_symbol", value = "n")

all_events_long$event_type <- unlist(lapply(all_events_long$gene_symbol, function(x)
    if (x %in% bg_events) {
        "Background"
    } else {
        if (x %in% pos_events) {
        "Positive"
    } else {
        "Negative"
    }
    }))

pos_events_long <- filter(all_events_long, event_type %in% c("Positive"))
neg_events_long <- filter(all_events_long, event_type %in% c("Negative"))

pos_events_long <- filter(all_events_long, event_type %in% c("Positive", "Background"))
neg_events_long <- filter(all_events_long, event_type %in% c("Negative", "Background"))
neg_events_long$n <- neg_events_long$n * -1

## zscores_by_pos <- spread(all_events_long, event_type, n) %>%
##     group_by(start) %>%
##         summarize(z_pos = (mean(Positive, na.rm = TRUE) - mean(Background, na.rm = TRUE)) / sd(Background, na.rm = TRUE),
##                   z_neg = (mean(Negative, na.rm = TRUE) - mean(Background, na.rm = TRUE)) / sd(Background, na.rm = TRUE))

t_by_pos <- spread(all_events_long, event_type, n) %>%
    group_by(start) %>%
        summarize(p_pos = t.test(Positive, Background, alternative = "g")$p.value,
                  p_neg = t.test(Negative, Background, alternative = "g")$p.value)

posPalette <- c("darkgrey", "red")
p_pos_motifs <- ggplot(pos_events_long) +
    geom_smooth(aes(y = n, x = start , color = event_type, fill = event_type), linetype = "longdash") +
        scale_colour_manual(values = posPalette) + scale_fill_manual(values = posPalette) +
            geom_point(data = filter(t_by_pos, p_pos < 0.01), aes(x = start, y = 0.6))

negPalette <- c("darkgrey", "blue")
p_neg_motifs <- ggplot(neg_events_long) +
    geom_smooth(aes(y = n, x = start , color = event_type, fill = event_type), linetype = "longdash") +
        scale_colour_manual(values = negPalette) + scale_fill_manual(values = negPalette) +
            geom_point(data = filter(t_by_pos, p_neg < 0.01), aes(x = start, y = -0.6))

print(p_neg_motifs)

tracks(Positive = p_pos_motifs, p_model, Negative = p_neg_motifs,
       xlim = c(-250,250),
       heights = c(1, 0.25, 1))



##----------------------------------
## n motifs averaged across positon
##----------------------------------
neg_event_motifs <- tbl_df(all_events[, neg_events])
pos_event_motifs <- tbl_df(all_events[, pos_events])
bg_event_motifs <- tbl_df(all_events[, bg_events])

averaged_events <- data_frame(start = all_events$start,
                              end = all_events$end,
                              chr = all_events$chr,
                              Negative = -rowMeans(neg_event_motifs),
                              Positive = rowMeans(pos_event_motifs),
                              Background = rowMeans(bg_event_motifs))

averaged_events <- gather(averaged_events, 4:6, key = "event_type", value = "n")

pos_events_long <- filter(averaged_events, event_type %in% c("Positive", "Background"))
neg_events_long <- filter(averaged_events, event_type %in% c("Negative", "Background"))

posPalette <- rev(c("darkgrey", "red"))
p_pos_motifs <- ggplot(pos_events_long) +
    geom_smooth(aes(y = n, x = start , color = event_type, fill = event_type), linetype = "longdash") +
        scale_colour_manual(values = posPalette) + scale_fill_manual(values = posPalette) +

negPalette <- rev(c("darkgrey", "blue"))
p_neg_motifs <- ggplot(neg_events_long) +
    geom_smooth(aes(y = n, x = start , color = event_type, fill = event_type), linetype = "longdash") +
        scale_colour_manual(values = negPalette) + scale_fill_manual(values = negPalette) +

tracks(Positive = p_pos_motifs, p_model, Negative = p_neg_motifs,
       xlim = c(-250,250),
       heights = c(1, 0.25, 1))

##----------------------------------
## using Gviz to plot motif distributions
##----------------------------------
library(Gviz)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

event <- "MBNL1"
event_info <- events_bed %>% filter(name == event)

chr <- unique(event_info$chr)
gen <- "hg19"

grtrack <- GeneRegionTrack(txdb, genome = gen,
                           chromosome = chr,
                           name = "UCSC known genes")

gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

lim <- c(min(event_info$start), max(event_info$end))

event_dat <- data_frame(coords = c(event_info[score == "E",]$start - as.numeric(names(motifs_per_event[[event]]$UI_Seq)),
                            event_info[score == "E",]$end + as.numeric(names(motifs_per_event[[event]]$DI_Seq))),
                        n_motifs = as.numeric(c(motifs_per_event[[event]]$UI_Seq, motifs_per_event[[event]]$DI_Seq)))

dtrack <- DataTrack(data = event_dat$n_motifs, start = event_dat$coords, width = 10,
                    chromosome = chr, genome = gen,
                    name = "# of MBNL motifs")

plotTracks(list(itrack, gtrack, grtrack, dtrack),
           from = lim[1], to = lim[2],
           type = c("p", "a"))
