getMaxPeakValue <- function(starts, stops, obj, ratio) {
  
  #get peaks for each base
  Apeaks <- getpeaks(obj@traceMatrix[,1])
  Cpeaks <- getpeaks(obj@traceMatrix[,2])
  Gpeaks <- getpeaks(obj@traceMatrix[,3])
  Tpeaks <- getpeaks(obj@traceMatrix[,4])
  #hack for last peak. Just uses distance preceding peak 
  #as distance after peak

 
  #Now get max peak value for each channel in each peak window. 
  #If no peak return 0 
  primary <- NULL
  secondary <- NULL
  tempPosMatrix <- matrix(nrow=length(starts), ncol=4)
  tempAmpMatrix <- matrix(nrow=length(starts), ncol=4)
  #obj <- setFrame(obj, starts)
  
  for(i in 1:length(starts)) {
    Apeak <- peakvalues(Apeaks, starts[i], stops[i])
    Cpeak <- peakvalues(Cpeaks, starts[i], stops[i])
    Gpeak <- peakvalues(Gpeaks, starts[i], stops[i])
    Tpeak <- peakvalues(Tpeaks, starts[i], stops[i])
    
    if(is.na(Apeak[2]) & 
       is.na(Cpeak[2]) & 
       is.na(Gpeak[2]) & 
       is.na(Tpeak[2])) next #rare case where no peak found 
    signals <- c(Apeak[1], Cpeak[1], Gpeak[1], Tpeak[1])
    tempAmpMatrix[i,] <- signals
    positions <- c(Apeak[2], Cpeak[2], Gpeak[2], Tpeak[2])
    tempPosMatrix[i,] <- positions
    signalratios <- signals/max(signals, na.rm=TRUE)
    Bases <- c("A", "C", "G", "T")
    Bases[signalratios < ratio] <- NA
    #sort by decreasing signal strength
    Bases <- Bases[order(signals, decreasing=TRUE)] 
    positions <- positions[order(signals, decreasing=TRUE)]
    
    if(length(Bases[!is.na(Bases)]) == 4 
       | length(Bases[!is.na(Bases)]) == 0) {
      primary <- c(primary, "N")
      secondary <- c(secondary, "N")
    }
    else if(length(Bases[!is.na(Bases)]) > 1) {
      primary <- c(primary, Bases[1]) 
      Bases2 <- Bases[2:4]
      secondary <- c(secondary, 
                     mergeIUPACLetters(paste(sort(Bases2[!is.na(Bases2)]), 
                                             collapse="")))
    }
    else {
      primary <- c(primary, Bases[1])
      secondary <- c(secondary, Bases[1])
    }
  }
  obj <- objectConsolidation(tempPosMatrix, tempAmpMatrix, DNAString, primary, secondary, obj)
  return(obj)
}

 objectConsolidation<- function(tempPosMatrix, tempAmpMatrix, DNAString, primary, secondary, obj) {
   
   obj@peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
   obj@peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
   obj@primarySeqID <- "sangerseq package primary basecalls"
   obj@primarySeq <- DNAString(paste(primary, collapse=""))
   obj@secondarySeqID <- "sangerseq package secondary basecalls"
   obj@secondarySeq <- DNAString(paste(secondary, collapse=""))
   return(obj)
 }
 

 
 #named vector in R (dictonary in python)