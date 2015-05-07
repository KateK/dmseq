##----------------------------------
## required libraries
##----------------------------------
library(data.table)
library(dplyr)
library(Biostrings)

##----------------------------------
## Sample Group
##----------------------------------
sample_group <- "tibialis"

##----------------------------------
## ALE
##----------------------------------
metadata_file <- "~/Projects/DMseq/data/ALE_metadata.txt"
results_file <- paste("~/Projects/DMseq/results/", sample_group, "/", sample_group, "_ALE_results.txt", sep = "")

##----------------------------------
## AFE
##----------------------------------
metadata_file <- "~/Projects/DMseq/data/AFE_metadata.txt"
results_file <- paste("~/Projects/DMseq/results/", sample_group, "/", sample_group, "_AFE_results.txt", sep = "")
                      

##----------------------------------
## Analysis Functions
##----------------------------------
analyzeUTRmotifs <- function(metadata_file, results_file)
{
    require(data.table)
    require(dplyr)
    require(Biostrings)

    metadata <- tbl_dt(fread(metadata_file))
    res <- tbl_dt(fread(results_file))

    sig_res <- res %>% filter(DM1_n_sig > 11,
                              abs(delta_psi_mean) >= 0.05)
    
    dysregulated_genes <- unique(sig_res$event_name)
    iso_seqs <- metadata %>% filter(gene_symbol %in% dysregulated_genes)

    mbnl_motifs <- paste(c("GCTT", "CGCT", "TGCT", "GCGC", "CCGC", "CTGC",
                           "GCTA", "ACGC", "CGCA", "AGCT", "TTGC", "CAGC"),
                         collapse = "|")
    celf_motifs <- paste(c("TGTT", "ATGT", "TTGT", "TGTC", "GTGT", "TGTA",
                           "GTTT", "TGTG", "GTCT", "TTTT"), collapse = "|")

    ## sliding window over sequence; counts motifs per chunk
    chunk_size <- 10
    foreach(i=1:2, .combine = list, .multicombine = TRUE) %do%
    {
        event <- iso_seqs[i,]
        starts <- seq(1, nchar(event$seq) - chunk_size, by = 1)
        n <- length(starts)
        motifs_per_chunk <- foreach(j=1:n, .combine = c, .multicombine = TRUE) %do%
        { 
            chunk <- substr(event$seq, starts[j], starts[j] + chunk_size -1)
            sum(str_count(mbnl_motifs, chunk))
            setNames(sum(str_count(chunk, mbnl_motifs)), starts[j])
        }
    }
    
    iso_seqs_obj <- DNAStringSet(iso_seqs$seq)
    names(iso_seqs_obj) <- iso_seqs$isoform

    ## a PWM for mbnl /  celf may be more appropriate... 
    mbnl_motifs <- DNAStringSet(c("GCTT", "CGCT", "TGCT", "GCGC", "CCGC", "CTGC",
                                  "GCTA", "ACGC", "CGCA", "AGCT", "TTGC", "CAGC"))
    celf_motifs <- DNAStringSet(c("TGTT", "ATGT", "TTGT", "TGTC", "GTGT", "TGTA",
                                  "GTTT", "TGTG", "GTCT", "TTTT"))

    lapply(mbnl_motifs, function(x) vmatchPattern(x, iso_seqs_obj))
    
}
