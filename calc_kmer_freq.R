library(stringi)
library(seqinr)
library(data.table)
library(gtools)

## Import utility functions
source("~/Analysis_Projects/DMseq/bin/kmer_tools.R")

setwd("~/Analysis_Projects/DMseq/")

utr_seqs <- fread("data/ALE_metadata.txt")
utr_seqs[, isoform_id := gsub("..$", "", isoform)]

test <- utr_seqs %>% head(100) %>% as.list
test$dist <- sequence_kmer_freq(test$sequence, 4)

sequence_kmer_freq <- function(sequence, k) {
    k_mers <- lapply(sequence, function(x) {
        seq_loop_size <- stri_length(x) - k + 1  
        kmers <- sapply(1:seq_loop_size, function(z) {
            y <- z + k - 1
            kmer <- substr(x = x, start = z, stop = y)
            return(kmer)
        })
        return(kmers)
    })
    uniq <- unique(unlist(k_mers))
    ind <- t(sapply(k_mers, function(x) {
        tabulate(match(x, uniq), length(uniq))
    }))
    colnames(ind) <- uniq    
    return(ind)
}
