library(stringi)
library(seqinr)
library(data.table)
library(gtools)

## Import utility functions
source("~/Projects/DMseq/bin/motif_utils.R")

setwd("~/Projects/DMseq/")

utr_seqs <- fread("data/ALE_metadata.txt")
utr_seqs[, isoform_id := gsub("..$", "", isoform)]

test <- utr_seqs %>% head(100) %>% as.list

test$dist <- lapply(lapply(test$sequence, getKmersFromSeq, 4), table)

test$dist

