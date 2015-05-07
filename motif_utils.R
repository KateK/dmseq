library(stringr)
library(gtools)

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
getKmersFromString <- function(string, k, unique=F) {
    seq_loop_size <- str_length(string) - k + 1  
    if (seq_loop_size >= 1) {
        kmers <- sapply(1:seq_loop_size, function(z) {
            y <- z + k - 1
            kmer <- substr(x = string, start = z, stop = y)
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
    if (str_length(query) > str_length(subject)) {
        return(numeric(0))
    } else {
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

##--------------------------------------------
## for each k-mer a in string
## for each k-mer a' differing from a by at most d mismatches
## if a' appears in each string from string with at most d mismatches
## output the number of times a' appears
##--------------------------------------------
kmerEnumerate <- function(stringCollection, k) {
    kmerized_strings <- lapply(stringCollection, getKmersFromString, k)
    enumerated_kmers <- lapply(kmerized_strings, table)    
    return(enumerated_kmers)
}

##--------------------------------------------
## for collection of kmers
## calculate the distance between each
## pair 
##
## output all clusters of kmers 
##--------------------------------------------
clusterKmers <- function(kmerCollection, maxDistance = 1) {

stop("ERROR: Not implemented")
    
}

##--------------------------------------------
## Find consensus for a set of kmers
##--------------------------------------------
findConsensus <- function(kmers, ...) {
    require(Biostrings)
    kmer_set <- DNAStringSet(kmers)
    consensusString(kmer_set, ...)    
}

##--------------------------------------------
## select random kmer from string
##--------------------------------------------
selectRandomKmer <- function(text, k){
    maxStart <- str_length(text) - k + 1
    randomKmerStart <- sample(1:maxStart, size=1)
    randomKmer <- substr(x=text, start=randomKmerStart, stop=randomKmerStart+k-1)
    return(randomKmer)
}

##--------------------------------------------
## convert set of motifs to counts matrix
##--------------------------------------------
motifsToCountMatrix <- function(motifs, pseudocount=F) {
    motifMatrix <- sapply(motifs, stringToCharVec)
    countMat <- apply(motifMatrix, 1, function(row)
        {
            freqA <- sum(row=="A")
            freqC <- sum(row=="C")
            freqG <- sum(row=="G")
            freqT <- sum(row=="T")
            return(c(freqA, freqC, freqG, freqT))
        })
    rownames(countMat) <- c("A","C","G","T")
    if (pseudocount) {
        countMat <- countMat + 1
    }
    return(countMat)
}

##--------------------------------------------
## convert set of motifs to prob matrix
##--------------------------------------------
motifsToProbMatrix <- function(motifs, pseudocount=F) {
    countMat <- motifsToCountMatrix(motifs, pseudocount)
    return(countMat / colSums(countMat))
}
