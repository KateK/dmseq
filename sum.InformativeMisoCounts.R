##' Calculate the total number of informative reads for a MISO event
##'
##' Informative reads in this context are those reads that consistent with a single isoform.
##'
##' For the two-isoform case in exon-centric analyses, the counts field general format is:
##' (1,0):X,(0,1):Y,(1,1):Z,(0,0):L
##' where X, Y, Z, L are integer counts corresponding to the number of reads in each of these categories. Class (1,0) are reads consistent with the first isoform in the annotation but not the second, class (0,1) are reads consistent with the second but not the first, class (1,1) are consistent with both isoforms, and reads in (0,0) are consistent with neither. 
##' @title sum.InformativeMisoCounts
##' @param counts_string character vector, each item pertaining to a MISO counts field for an event
##' @return numeric vector
##' @author Adam Struck
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
