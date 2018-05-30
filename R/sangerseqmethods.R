#Implements the show generic 
setMethod("show", "sangerseq",
  function(object) {
    cat("Number of datapoints: ")
    cat(nrow(object@traceMatrix))
    cat("\n")
    cat("Number of basecalls: ")
    cat(length(object@primarySeq))
    cat("\n\n")
    cat("Primary Basecalls: ")
    cat(toString(object@primarySeq))
    cat("\n\n")
    cat("Secondary Basecalls: ")
    cat(toString(object@secondarySeq))
    cat("\n")
  }
)

#Constructors
#' @rdname sangerseq-class
#' @aliases sangerseq,abif-method
setMethod("sangerseq", "abif", 
  function(obj) {
    res <- new("sangerseq")
    tracematrix <- matrix(c(obj@data$DATA.9, 
                            obj@data$DATA.10, 
                            obj@data$DATA.11, 
                            obj@data$DATA.12), 
                          ncol=4)
    orderedmatrix <- cbind(tracematrix[,regexpr("A", obj@data$FWO_.1)[1]],
                           tracematrix[,regexpr("C", obj@data$FWO_.1)[1]],
                           tracematrix[,regexpr("G", obj@data$FWO_.1)[1]],
                           tracematrix[,regexpr("T", obj@data$FWO_.1)[1]]
    )
    

    #check for valid chars
    basecalls1 <- strsplit(obj@data$PBAS.2, "")[[1]]
    basecalls1 <- paste0(basecalls1[basecalls1 %in% DNA_ALPHABET], 
                         collapse = "")
    if (nchar(basecalls1) != nchar(obj@data$PBAS.2)) {
      warning("Invalid characters removed from primary basecalls. This may result
                in basecalls being shifted. Please check chromatogram.")
    }
    #Appears normal to have them not match
    #if (nchar(basecalls1) != length(obj@data$PLOC.2)) {
    #  warning("Number of primary basecalls does not match the number of peaks. Please
    #          check chromatogram.")
    #}

    basecalls1 <- DNAString(substr(basecalls1,1,length(obj@data$PLOC.2)))
    basecallpositions1 <- obj@data$PLOC.2 + 1
    if(!is.null(obj@data$P2BA.1)) {
      basecalls2 <- strsplit(obj@data$P2BA.1, "")[[1]]
      basecalls2 <- paste0(basecalls2[basecalls2 %in% DNA_ALPHABET], 
                           collapse = "")
      if (nchar(basecalls2) != nchar(obj@data$P2BA.1)) {
        warning("Invalid characters removed from secondary basecalls. This may 
                result in basecalls being shifted. Please check chromatogram.")
      }
      #Appears normal to have them not match
      #if (nchar(basecalls2) != length(obj@data$PLOC.2)) {
      #  warning("Number of secondary basecalls does not match the number of peaks. Please
      #        check chromatogram.")
      #}
      basecalls2 <- DNAString(substr(basecalls2,1,length(obj@data$PLOC.2)))
      basecallpositions2 <-obj@data$PLOC.2 + 1
    } else {
      basecalls2 <- DNAString("")
      basecallpositions2 <- NA
    }
    if(!is.null(obj@data$P1AM.1)) {
      peakamps1 <- obj@data$P1AM.1
    } else {
      peakamps1 <- NA
    }
    if(!is.null(obj@data$P2AM.1)) {
      peakamps2 <- obj@data$P2AM.1
    } else {
      peakamps2 <- NA
    }
    res@primarySeqID <- "From ab1 file"
    res@primarySeq <- basecalls1
    res@secondarySeqID <- "From ab1 file"
    res@secondarySeq <- basecalls2
    res@traceMatrix <- orderedmatrix
    res@peakPosMatrix <- cbind(basecallpositions1, basecallpositions2, 
                               NA, NA, deparse.level=0)
    res@peakAmpMatrix <- cbind(peakamps1, peakamps2, NA, NA, 
                               deparse.level=0)
    return(res)
  }
)


#' @rdname sangerseq-class
#' @aliases sangerseq,scf-method
setMethod("sangerseq", "scf", 
  function(obj) {
    res <- new("sangerseq")
    res@primarySeqID <- "From scf file"
    res@primarySeq <- DNAString(obj@basecalls)
    res@secondarySeqID <- "From scf file"
    res@secondarySeq <- DNAString("")
    res@traceMatrix <- obj@sample_points
    res@peakPosMatrix <- cbind(obj@basecall_positions, NA, NA, NA, 
                               deparse.level=0)
    #res@peakAmpMatrix <- nothing assigned because data not in file
    return(res)
  }
)

#' @rdname makeBaseCalls
setMethod("makeBaseCalls", "sangerseq",
  function(obj, ratio=.33) {
    
    #get peaks for each base
    Apeaks <- getpeaks(obj@traceMatrix[,1])
    Cpeaks <- getpeaks(obj@traceMatrix[,2])
    Gpeaks <- getpeaks(obj@traceMatrix[,3])
    Tpeaks <- getpeaks(obj@traceMatrix[,4])
    
    #get window around primary basecall peaks
    primarypeaks <- obj@peakPosMatrix[,1]
    diffs <- diff(c(0,primarypeaks))
    starts <- primarypeaks - 0.5*diffs
    stops <- c(primarypeaks[1:(length(primarypeaks)-1)] + 
                 0.5*diffs[2:length(diffs)], 
               primarypeaks[length(diffs)] + 0.5*diffs[length(diffs)]
    ) 
    #hack for last peak. Just uses distance preceding peak 
    #as distance after peak
    
    #Now get max peak value for each channel in each peak window. 
    #If no peak return 0  
    primary <- NULL
    secondary <- NULL
    tempPosMatrix <- matrix(nrow=length(starts), ncol=4)
    tempAmpMatrix <- matrix(nrow=length(starts), ncol=4)
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
    obj@peakPosMatrix <- tempPosMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    obj@peakAmpMatrix <- tempAmpMatrix[rowSums(!is.na(tempPosMatrix)) > 0,]
    obj@primarySeqID <- "sangerseq package primary basecalls"
    obj@primarySeq <- DNAString(paste(primary, collapse=""))
    obj@secondarySeqID <- "sangerseq package secondary basecalls"
    obj@secondarySeq <- DNAString(paste(secondary, collapse=""))
    
    return(obj)
  }
)

#' @rdname chromatogram
setMethod("chromatogram", "sangerseq", 
  function(obj, trim5=0, trim3=0, 
           showcalls=c("primary", "secondary", "both", "none"), 
           width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
           filename=NULL, showtrim=FALSE, showhets=TRUE) {
    originalpar <- par(no.readonly=TRUE)
    showcalls <- showcalls[1]
    traces <- obj@traceMatrix
    basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
    basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
    averagePosition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
    basecalls1 <- basecalls1[1:length(averagePosition)] 
    basecalls2 <- basecalls2[1:length(averagePosition)] 
    
    if(showtrim == FALSE) {
      trimResult <- removeTrim(basecalls1, basecalls2, trim5, trim3)
      basecalls1 <- trimResult[["basecalls1"]]
      basecalls2 <- trimResult[["basecalls2"]]
      averagePosition <- averagePosition[(1 + trim5):(length(averagePosition) - trim3)] 
    }
    
    settings <- makeSettings(basecalls1, basecalls2, trim5, trim3, 
                             averagePosition, traces, ylim, width)
    
    traces <- settings[["traces"]]
    basecalls1 <- settings[["basecalls1"]] 
    basecalls2 <- settings[["basecalls2"]]
    averagePosition <- settings[["averagePosition"]]
    ylims <- settings[["ylims"]]
    startTrims <- settings[["startTrims"]]
    endTrims <- settings[["endTrims"]] 
    startHets <- settings[["startHets"]]
    endHets <- settings[["endHets"]]
    colorVector1 <- settings[["colorVector1"]] 
    colorVector2 <- settings[["colorVector2"]]
    traceWidth <- settings[["traceWidth"]]
    breaks <- settings[["breaks"]]
    numPlots <- settings[["numPlots"]]
    
    if(!is.null(filename)) {
      pdf(filename, width=8.5, height=height*numPlots)
    }
    par(mar=c(2,2,2,1), mfrow=c(numPlots, 1))
    basecallwarning1 = 0
    basecallwarning2 = 0
    j = 1
    
    for(i in breaks) {
      
      range <- averagePosition >= i & averagePosition < (i+traceWidth)
      startHet <- startHets[range] - traceWidth*(j-1)
      startHet[startHet < 0] <- 0
      endHet <- endHets[range] - traceWidth*(j-1)
      endHet[endHet > traceWidth] <- traceWidth
      lab1 <- basecalls1[range]
      lab2 <- basecalls2[range]
      pos <- averagePosition[range] - traceWidth*(j-1)
      colors1 <- colorVector1[range]
      colors2 <- colorVector2[range]
      startTrim <- startTrims[range] - traceWidth*(j-1)
      endTrim <- endTrims[range] - traceWidth*(j-1)
      plotRange <- i:min(i+traceWidth, nrow(traces))
      plot(traces[plotRange,1], type='n', ylim=ylims, ylab="", xaxt="n", 
           bty="n", xlab="", yaxt="n", , xlim=c(1,traceWidth))
      if (showhets==TRUE) {
        rect(startHet, 0, endHet, ylims[2], col='#D5E3F7', border='#D5E3F7')
      }
      if (showtrim==TRUE) {
        rect(startTrim, 0, endTrim, ylims[2], col='red', border='transparent', 
             density=15)
      }
      lines(traces[plotRange,1], col="green")
      lines(traces[plotRange,2], col="blue")
      lines(traces[plotRange,3], col="black")
      lines(traces[plotRange,4], col="red")
      mtext(as.character(which(range)[1]), side=2, line=0, cex=cex.mtext)
      
      for(k in 1:length(lab1)) {
        if (showcalls=="primary" | showcalls=="both") {
          if (is.na(basecalls1[1]) & basecallwarning1==0) {
            warning("Primary basecalls missing")
            basecallwarning1 = 1
          } 
          else if (length(lab1) > 0) {   
            axis(side=3, at=pos[k], labels=lab1[k], col.axis=colors1[k], 
                 family="mono", cex=cex.base, line=ifelse(showcalls=="both", 0, 
                                                          -1), tick=FALSE)
          }
        }
        if (showcalls=="secondary" | showcalls=="both") {
          if (is.na(basecalls2[1]) & basecallwarning2 == 0) {
            warning("Secondary basecalls missing")
            basecallwarning2 = 1
          } 
          else if (length(lab2) > 0) { 
            axis(side=3, at=pos[k], labels=lab2[k], col.axis=colors2[k], 
                 family="mono", cex=cex.base, line=-1, tick=FALSE)
          }
        }
      }
      j = j + 1
    }
    
    if(!is.null(filename)) {
      dev.off()
      cat(paste("Chromatogram saved to", filename, 
                "in the current working directory"))
    }
    else par(originalpar)
  }
)

#' @rdname setAllelePhase
setMethod("setAllelePhase", "sangerseq",
  function(obj, refseq, trim5 = 0, trim3 = 0, asIs = FALSE) {
    
    if (!asIs) {
      obj <- makeBaseCalls(obj)
    }
    
    refseq <- DNAString(refseq)
    basecalls <- basecalldf(obj)
    
    refstrandresult <- determinerefstrand(refseq, basecalls, trim5, trim3)
    refstrand <- refstrandresult[["refstrand"]]
    pairwisealigned <- refstrandresult[["pa"]]

    #This picks the longest continuous alignment range
    if(pairwisealigned@pattern@range@width < 10) {
      stop("Seed length not long enough. Must provide at least 10 bases of good 
           matching sequence.")
    }
    
    basecalls <- accountforstartoffset(refseq, pairwisealigned, basecalls, trim5)
    basecalls <- updateprimarysecondary(basecalls)
    
    obj <- updatesangerseqobj(obj, refstrand, basecalls)
    return(obj)
  }
)

#accessors
#' @rdname sangerseqAccessors
#' @aliases primarySeq

setMethod("primarySeq", "sangerseq",
  function(obj, string=FALSE) {
    if(string==TRUE) toString(obj@primarySeq)
    else obj@primarySeq
  }
)

#' @rdname sangerseqAccessors
#' @aliases secondarySeq

setMethod("secondarySeq", "sangerseq",
  function(obj, string=FALSE) {
    if(string==TRUE) toString(obj@secondarySeq)
    else obj@secondarySeq
  }
)

#' @rdname sangerseqAccessors
#' @aliases traceMatrix

setMethod("traceMatrix", "sangerseq", function(obj) obj@traceMatrix)

#' @rdname sangerseqAccessors
#' @aliases peakPosMatrix

setMethod("peakPosMatrix", "sangerseq", function(obj) obj@peakPosMatrix)

#' @rdname sangerseqAccessors
#' @aliases peakAmpMatrix

setMethod("peakAmpMatrix", "sangerseq", function(obj) obj@peakAmpMatrix)

#' @rdname sangerseqAccessors
#' @aliases primarySeqID

setMethod("primarySeqID", "sangerseq", function(obj) obj@primarySeqID)

#' @rdname sangerseqAccessors
#' @aliases secondarySeqID

setMethod("secondarySeqID", "sangerseq", function(obj) obj@secondarySeqID)

#setters
#' @rdname sangerseqAccessors
#' @aliases primarySeq<-

setMethod("primarySeq<-", "sangerseq",
  function(obj, value) {
    if(class(value)!="DNAString") value <- DNAString(value)
    obj@primarySeq <- value 
    obj
  })

#' @rdname sangerseqAccessors
#' @aliases secondarySeq<-

setMethod("secondarySeq<-", "sangerseq",
  function(obj, value) {
    if(class(value)!="DNAString") value <- DNAString(value)
    obj@secondarySeq <- value 
    obj
  })

#' @rdname sangerseqAccessors
#' @aliases traceMatrix<-

setMethod("traceMatrix<-", "sangerseq",
  function(obj, value) {
    obj@traceMatrix <- value 
    obj
  })

#' @rdname sangerseqAccessors
#' @aliases peakPosMatrix<-

setMethod("peakPosMatrix<-", "sangerseq",
  function(obj, value) {
    obj@peakPosMatrix <- value 
    obj
  })

#' @rdname sangerseqAccessors
#' @aliases peakAmpMatrix<-

setMethod("peakAmpMatrix<-", "sangerseq",
  function(obj, value) {
    obj@peakAmpMatrix <- value 
    obj
  })

#' @rdname sangerseqAccessors
#' @aliases primarySeqID<-

setMethod("primarySeqID<-", "sangerseq",
  function(obj, value) {
    obj@primarySeqID <- value 
    obj
  })

#' @rdname sangerseqAccessors
#' @aliases secondarySeqID<-

setMethod("secondarySeqID<-", "sangerseq",
  function(obj, value) {
    obj@secondarySeqID <- value 
    obj
  })