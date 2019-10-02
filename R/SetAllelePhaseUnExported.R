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

updatePrimarySecondary <- function(obj) {
  
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


determineRefStrand <- function(refSeq, basecalls, trim5, trim3) {
  
  seedSeq <- paste(basecalls$consensus[(trim5+1):(nrow(basecalls)-trim3)], 
                   collapse="")
  
  pa <- pairwiseAlignment(seedSeq, refSeq, type="local", 
                          gapOpening=-200, gapExtension=-10)
  paRC <- pairwiseAlignment(seedSeq, reverseComplement(DNAString(refSeq)), 
                            type="local", gapOpening=-200, gapExtension=-10)
  refStrand <- "Allele 2"
  #refStrand <- "Reference"
  if(score(paRC) > score(pa)) {
    refSeq <- toString(reverseComplement(refSeq))
    pa <- paRC
    refStrand <- "Reference (revcomp)"
  }
  
  return(list(refStrand = refStrand, pa = pa))
}


accountForStartOffset <- function(refSeq, pa, basecalls, trim5) {
  
  refVector <- strsplit(toString(refSeq), "")[[1]]
  seqStart <- pa@pattern@range@start 
  refStart <- pa@subject@range@start
  startOffset <- refStart - seqStart - trim5
  if (startOffset < 0)  {
    refVector <- c(rep("N", abs(startOffset)), refVector)
    startOffset <- 0
  }
  end <- nrow(basecalls) + startOffset
  if (end > length(refVector)) {
    refVector <- c(refVector, rep("N", end-length(refVector)))
  }
  if ((refStart - seqStart - trim5) < 0 | end > length(refVector)) {
    warning("Reference sequence does not encompass sequencing results. 
            Ambiguous bases will be attributed to both alleles outside the region 
            covered by the reference sequence.\n")
  }  
  
  basecalls$ref <- refVector[(startOffset + 1):end]
  
  return(basecalls)
}


updateSangerSeqObj <- function(obj, refStrand, basecalls) {
  
  #Update Sangerseq Obj
  obj@primarySeqID <- refStrand
  primarySeq <- paste0(mergeIUPACLetters(basecalls$newprimary), collapse="")
  primarySeq(obj) <- DNAString(primarySeq)
  obj@secondarySeqID <- paste("NonReference")
  secondarySeq <- paste0(mergeIUPACLetters(basecalls$newsecondary), 
                         collapse="")
  secondarySeq(obj) <- DNAString(secondarySeq)
  return(obj)
}