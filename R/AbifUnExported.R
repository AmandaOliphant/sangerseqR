makeOrderedMatrix <- function(obj) {
  
  traceMatrix <- matrix(c(obj@data$DATA.9, 
                          obj@data$DATA.10, 
                          obj@data$DATA.11, 
                          obj@data$DATA.12), 
                        ncol=4)
  orderedMatrix <- cbind(traceMatrix[,regexpr("A", obj@data$FWO_.1)[1]],
                         traceMatrix[,regexpr("C", obj@data$FWO_.1)[1]],
                         traceMatrix[,regexpr("G", obj@data$FWO_.1)[1]],
                         traceMatrix[,regexpr("T", obj@data$FWO_.1)[1]])
  
  return(orderedMatrix)
  
}

getBasecalls <- function(basecallData, peakLocations) {
  
  #check for valid chars
  basecalls <- strsplit(basecallData, "")[[1]]
  basecalls <- paste0(basecalls[basecalls %in% DNA_ALPHABET], 
                       collapse = "")
  if (nchar(basecalls) != nchar(basecallData)) {
    warning("Invalid characters removed from primary basecalls. This may result
            in basecalls being shifted. Please check chromatogram.")
  }
  #Appears normal to have them not match
  #if (nchar(basecalls1) != length(obj@data$PLOC.2)) {
  #  warning("Number of primary basecalls does not match the number of peaks. Please
  #          check chromatogram.")
  #}
  
  basecalls <- DNAString(substr(basecalls, 1, length(peakLocations)))
  
  return(basecalls)
  
}

makePeakAmpMatrix <- function(obj) {
  
  if(!is.null(obj@data$P1AM.1)) {
    peakamps1 <- obj@data$P1AM.1
  } 
  else {
    peakamps1 <- NA
  }
  
  if(!is.null(obj@data$P2AM.1)) {
    peakamps2 <- obj@data$P2AM.1
  } 
  else {
    peakamps2 <- NA
  }
  
  newPeakAmpMatrix <- cbind(peakamps1, peakamps2, NA, NA, 
                             deparse.level=0)
  
  return(newPeakAmpMatrix)
}