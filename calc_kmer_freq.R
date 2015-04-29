library(stringi)
library(seqinr)

setwd("~/Projects/DMseq/")

se_seqs <- read.fasta(file = "data/SE.cis500.fa", seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE)

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
