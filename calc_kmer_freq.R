library(data.table)
library(seqinr)
library(gtools)
library(dplyr)

## Import utility functions
source("~/Projects/DMseq/bin/motif_utils.R")

setwd("~/Projects/DMseq/")

utr_seqs <- fread("data/ALE_metadata.txt")
utr_seqs[, isoform_id := gsub("..$", "", isoform)]

test <- utr_seqs %>% head(100) %>% as.list

test$dist <- kmerEnumerate(test$sequence, 6)
test$dist
