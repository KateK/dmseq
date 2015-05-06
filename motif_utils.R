library("stringr")
library("gtools")

##--------------------------------------------
## Get all possible DNA Kmers 
##--------------------------------------------
getAllPossibleKmers <- function(len) {
    allDnaVector <- permutations(n=4, r=len,
                                 v=c("A","T","G","C"),
                                 repeats.allowed=TRUE)
    allDna <- apply(allDnaVector, 1, paste0, collapse="")
    return(allDna)
}

##--------------------------------------------
## given some sequence, find all possible kmers 
##--------------------------------------------
getKmersFromSeq <- function(sequence, k, unique=F) {
    seq_loop_size <- str_length(sequence) - k + 1  
    if (seq_loop_size >= 1) {
        kmers <- sapply(1:seq_loop_size, function(z) {
            y <- z + k - 1
            kmer <- substr(x = sequence, start = z, stop = y)
            return(kmer)                
        })
    } else {
        kmers <- NA
    }            
    if (unique) {
        return(unique(kmers))
    } else {
        return(kmers)
    }
}

##--------------------------------------------
## Given a kmer, Find all kmers differing by
## at most 'd' mutations
##--------------------------------------------
getApproximateKmers <- function(kmer, d) {
    if (d == 0) {
        return(kmer)
    }
    kmerVec <- strsplit(kmer, "")[[1]]
    changedPos <- combn(1:length(kmerVec), d)
    changedLetters <- permutations(4, d, c("A","T","G","C"), repeats.allowed=T)   
    allPositionsApproxKmers <- sapply(1:ncol(changedPos), function(changePosNo) {
        onePositionApproxKmers <- sapply(1:nrow(changedLetters), function(changeLettersNo) {
            approxKmerVec <- kmerVec
            approxKmerVec[changedPos[,changePosNo]] <- changedLetters[changeLettersNo,]
            approxKmer <- paste0(approxKmerVec, collapse="")
            return(approxKmer)
        })        
        return(onePositionApproxKmers)
    })
    return(unique(as.vector(allPositionsApproxKmers)))
}

##--------------------------------------------
## convert a string to C-like char vector 
##--------------------------------------------
stringToCharVec <- function(string) {
    if(is.null(string) || is.na(string))
        return(c())

    charVec <- strsplit(string, split="")[[1]]
    return(charVec)
}

##--------------------------------------------
## calculate hamming distance between
## 2 strings
##--------------------------------------------
hammingDistance <- function(query, subject) {
    ## Only allow one query & subject
    stopifnot(length(subject) == 1 & length(query) == 1)
    ## make sure query is shorter than subject
    stopifnot(str_length(query) < str_length(subject))
    ## Define starts and ends based on query length
    start <- 1:(str_length(subject) - str_length(query) + 1)
    end <- str_length(query):str_length(subject)
    ## calculate distance
    ch1 <- stringToCharVec(query)
    ch2 <- substring(subject, start, end)    
    ch2 <- lapply(ch2, stringToCharVec)
    res <- unlist(lapply(ch2, function(x) sum(ch1 != x)))
    return(setNames(res, start))
}

##--------------------------------------------
## Given a dna (string), a segment (pattern) 
## Find all the positions in the dna  
## where the segment matches the dna  
## (at most d mutations is acceptable)
##--------------------------------------------
findMatches <- function(pattern, string, d=0, indexBase=1) {
    ## Hamming distance is number of mismatches
    potential_matches <- hammingDistance(pattern, string) <= d
    matchIndexes <- which(potential_matches==TRUE)
    ## Offset coordinates
    if(indexBase != 1) {
        matchIndexes <- matchIndexes - 1 + indexBase
    }    
    return(as.numeric(matchIndexes))
}

## MOTIFENUMERATION(Dna, k, d)
## for each k-mer a in Dna
## for each k-mer a' differing from a by at most d mutations
## if a' appears in each string from Dna with at most d mutations
## output the number of times a' appears
motifEnumerate <- function(dnaCollection, k, d) {
    allKmers <- lapply(dnaCollection, getKmersFromSeq, k, unique = TRUE)
    allKmers <- unique(unlist(allKmers))

    allApproximateKmers <- lapply(allKmers, getApproximateKmers, d)
    allApproximateKmers <- unique(unlist(allApproximateKmers))
   
    hasSegmentInDna <- function(dna, segment, d=0) {
        matchIndexes <- findMatches(string=dna, pattern=segment, d=d)
        return(length(matchIndexes))
    }

    hasSegmentInAllDna <- function(dnaCollection, segment, d=0) {
        hasInDnaCollection <- sapply(dnaCollection, hasSegmentInDna, segment, d)
        return(hasInDnaCollection)
    }

    isMotif <- sapply(allApproximateKmers, function(kmer)
        {
            return(hasSegmentInAllDna(dnaCollection, kmer, d))
        })
    
    motifs <- allApproximateKmers[isMotif]
    return(motifs)
}

# debug(motifEnumerate)
