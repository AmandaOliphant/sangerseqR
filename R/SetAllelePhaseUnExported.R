#Create data frame of basecalls
basecalldf <- function(obj) {
  primary <- strsplit(toupper(primarySeq(obj, string=TRUE)), "")[[1]]
  secondary <- strsplit(toupper(secondarySeq(obj, string=TRUE)), "")[[1]]
  basecalls <- data.frame(primary=primary, 
                          secondary=secondary, 
                          stringsAsFactors=FALSE)
  basecalls$primary <- unname(IUPAC_CODE_MAP[basecalls$primary])
  basecalls$secondary <- unname(IUPAC_CODE_MAP[basecalls$secondary])
  basecalls$consensus <- basecalls$primary
  basecalls$consensus[basecalls$primary != basecalls$secondary 
                      | nchar(basecalls$consensus) > 1] <- "N"
  basecalls$possibilities <- basecalls$consensus
  basecalls$possibilities[basecalls$possibilities == "N"] <- 
    paste0(basecalls$primary[basecalls$possibilities == "N"], 
           basecalls$secondary[basecalls$possibilities == "N"])
  return(basecalls)
}

updateprimarysecondary <- function(obj) {
  
  obj$noref <- apply(obj, 1, function(x) gsub(paste("[", x[5], "]"), "", x[4]))
  #if removing reference did not change nchar, then we still don't know so put
  #sequenced base(s) for both
  refNotOption <- nchar(obj$noref) == nchar(obj$possibilities)
  obj$newprimary[refNotOption] <- 
    obj$possibilities[refNotOption]
  obj$newsecondary[refNotOption] <- 
    obj$possibilities[refNotOption]
  #if not equal then ref goes in primary and noref (i.e. other(s)) goes in
  #secondary
  refNotOnlyOption <- nchar(obj$noref) != nchar(obj$possibilities)
  obj$newprimary[refNotOnlyOption] <- 
    obj$ref[refNotOnlyOption]
  obj$newsecondary[refNotOnlyOption] <- 
    obj$noref[refNotOnlyOption]
  #finally, if ref was only base, then secondary now has empty string. Replace
  #with ref
  refOnlyOption <- nchar(obj$newsecondary) == 0
  obj$newsecondary[refOnlyOption] <- obj$ref[refOnlyOption]
  
  return(obj)
}


determinerefstrand <- function(refseq, basecalls, trim5, trim3) {
  
  seedseq <- paste(basecalls$consensus[(trim5+1):(nrow(basecalls)-trim3)], 
                   collapse="")
  
  pa <- pairwiseAlignment(seedseq, refseq, type="local", 
                          gapOpening=-200, gapExtension=-10)
  paRC <- pairwiseAlignment(seedseq, reverseComplement(DNAString(refseq)), 
                            type="local", gapOpening=-200, gapExtension=-10)
  refstrand <- "Reference"
  if(score(paRC) > score(pa)) {
    refseq <- toString(reverseComplement(refseq))
    pa <- paRC
    refstrand <- "Reference (revcomp)"
  }
  
  return(list(refstrand = refstrand, pa = pa))
}


accountforstartoffset <- function(refseq, pa, basecalls, trim5) {
  
  refvector <- strsplit(toString(refseq), "")[[1]]
  seqstart <- pa@pattern@range@start 
  refstart <- pa@subject@range@start
  startoffset <- refstart - seqstart - trim5
  if (startoffset < 0)  {
    refvector <- c(rep("N", abs(startoffset)), refvector)
    startoffset <- 0
  }
  end <- nrow(basecalls) + startoffset
  if (end > length(refvector)) {
    refvector <- c(refvector, rep("N", end-length(refvector)))
  }
  if ((refstart - seqstart - trim5) < 0 | end > length(refvector)) {
    warning("Reference sequence does not encompass sequencing results. 
            Ambiguous bases will be attributed to both alleles outside the region 
            covered by the reference sequence.\n")
  }  
  
  basecalls$ref <- refvector[(startoffset + 1):end]
  
  return(basecalls)
}


updatesangerseqobj <- function(obj, refstrand, basecalls) {
  
  #Update Sangerseq Obj
  obj@primarySeqID <- refstrand
  primaryseq <- paste0(mergeIUPACLetters(basecalls$newprimary), collapse="")
  primarySeq(obj) <- DNAString(primaryseq)
  obj@secondarySeqID <- paste("NonReference")
  secondaryseq <- paste0(mergeIUPACLetters(basecalls$newsecondary), 
                         collapse="")
  secondarySeq(obj) <- DNAString(secondaryseq)
  return(obj)
}