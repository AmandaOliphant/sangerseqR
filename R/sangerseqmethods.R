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

    orderedMatrix <- makeOrderedMatrix(obj)
    
    basecalls1 <- getBasecalls(obj@data$PBAS.2, obj@data$PLOC.2)
    basecallpositions1 <- obj@data$PLOC.2 + 1
    
    if(!is.null(obj@data$P2BA.1)) {
     basecalls2 = getBasecalls(obj@data$P2BA.1, obj@data$PLOC.2)
      basecallpositions2 <-obj@data$PLOC.2 + 1
    } 
    else {
      basecalls2 <- DNAString("")
      basecallpositions2 <- NA
    }
    
    res@primarySeqID <- "From ab1 file"
    res@primarySeq <- basecalls1
    res@secondarySeqID <- "From ab1 file"
    res@secondarySeq <- basecalls2
    res@traceMatrix <- orderedMatrix
    res@peakAmpMatrix <- makePeakAmpMatrix(obj)
    res@peakPosMatrix <- cbind(basecallpositions1, basecallpositions2, 
                               NA, NA, deparse.level=0)
    
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
    
    #get window around primary basecall peaks
    primarypeaks <- obj@peakPosMatrix[,1]
    diffs <- diff(c(0,primarypeaks))
    starts <- primarypeaks - 0.5*diffs
    stops <- c(primarypeaks[1:(length(primarypeaks)-1)] + 
                 0.5*diffs[2:length(diffs)], 
               primarypeaks[length(diffs)] + 0.5*diffs[length(diffs)]
    ) 
     
    return(getMaxPeakValue(starts, stops, obj, ratio))
    
  }
)

setMethod("chromatogram", "character",
          function(obj, trim5=0, trim3=0, 
                   showcalls=c("primary", "secondary", "both", "none"), 
                   width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
                   filename=NULL, showtrim=FALSE, showhets=TRUE) {
            newSangerSeq <- readsangerseq(obj)
            chromatogram(newSangerSeq)
          }
)

setMethod("chromatogram", "abif",
  function(obj, trim5=0, trim3=0, 
          showcalls=c("primary", "secondary", "both", "none"), 
          width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
          filename=NULL, showtrim=FALSE, showhets=TRUE) {
    newSangerSeq <- sangerseq(obj)
    chromatogram(newSangerSeq)
  }
)

setMethod("chromatogram", "scf",
          function(obj, trim5=0, trim3=0, 
                   showcalls=c("primary", "secondary", "both", "none"), 
                   width=100, height=2, cex.mtext=1, cex.base=1, ylim=3, 
                   filename=NULL, showtrim=FALSE, showhets=TRUE) {
            newSangerSeq <- sangerseq(obj)
            chromatogram(newSangerSeq)
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
            aveposition <- rowMeans(obj@peakPosMatrix, na.rm=TRUE)
            basecalls1 <- basecalls1[1:length(aveposition)] 
            basecalls2 <- basecalls2[1:length(aveposition)] 
            if(showtrim == FALSE) {
              if(trim5+trim3 > length(basecalls1)) basecalls1 <- ""
              else basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
              if(trim5+trim3 > length(basecalls2)) basecalls2 <- ""
              else basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
              aveposition <- aveposition[(1 + trim5):(length(aveposition) - trim3)] 
            }
            indexes <- 1:length(basecalls1)
            trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all 
            #false if not trimmed
            if (!is.null(trim3)) {
              traces <- traces[1:(min(max(aveposition, na.rm=TRUE) + 10, 
                                      nrow(traces))), ]
            }
            if (!is.null(trim5)) {
              offset <- max(c(1, aveposition[1] - 10))
              traces <- traces[offset:nrow(traces),]
              aveposition <- aveposition - (offset-1)
            }
            maxsignal <- apply(traces, 1, max)
            ylims <- c(0, quantile(maxsignal, .75)+ylim*IQR(maxsignal))           
            p <- c(0, aveposition, nrow(traces))
            midp <- diff(p)/2
            starts <- aveposition - midp[1:(length(midp)-1)]
            starthets <- starts
            starthets[basecalls1 == basecalls2] <- NA
            ends <- aveposition + midp[2:(length(midp))]
            endhets <- ends
            endhets[basecalls1 == basecalls2] <- NA
            starttrims <- starts
            starttrims[!trimmed] <- NA
            endtrims <- ends
            endtrims[!trimmed] <- NA
            
            colortranslate <- c(A="green", C="blue", G="black", T="red")
            colorvector1 <- unname(colortranslate[basecalls1])
            colorvector1[is.na(colorvector1)] <- "purple"
            colorvector2 <- unname(colortranslate[basecalls2])
            colorvector2[is.na(colorvector2)] <- "purple"
            
            valuesperbase <- nrow(traces)/length(basecalls1)
            tracewidth <- width*valuesperbase
            breaks <- seq(1,nrow(traces), by=tracewidth)
            numplots <- length(breaks)
            if(!is.null(filename)) pdf(filename, width=8.5, height=height*numplots) 
            par(mar=c(2,2,2,1), mfrow=c(numplots, 1))
            basecallwarning1 = 0
            basecallwarning2 = 0
            j = 1
            
            for(i in breaks) {
              range <- aveposition >= i & aveposition < (i+tracewidth)
              starthet <- starthets[range] - tracewidth*(j-1)
              starthet[starthet < 0] <- 0
              endhet <- endhets[range] - tracewidth*(j-1)
              endhet[endhet > tracewidth] <- tracewidth
              lab1 <- basecalls1[range]
              lab2 <- basecalls2[range]
              pos <- aveposition[range] - tracewidth*(j-1)
              colors1 <- colorvector1[range]
              colors2 <- colorvector2[range]
              starttrim <- starttrims[range] - tracewidth*(j-1)
              endtrim <- endtrims[range] - tracewidth*(j-1)
              plotrange <- i:min(i+tracewidth, nrow(traces))
              plot(traces[plotrange,1], type='n', ylim=ylims, ylab="", xaxt="n", 
                   bty="n", xlab="", yaxt="n", , xlim=c(1,tracewidth))
              if (showhets==TRUE) {
                rect(starthet, 0, endhet, ylims[2], col='#D5E3F7', border='#D5E3F7')
              }
              if (showtrim==TRUE) {
                rect(starttrim, 0, endtrim, ylims[2], col='red', border='transparent', 
                     density=15)
              }
              lines(traces[plotrange,1], col="green")
              lines(traces[plotrange,2], col="blue")
              lines(traces[plotrange,3], col="black")
              lines(traces[plotrange,4], col="red")
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
  function(obj, refSeq, trim5 = 0, trim3 = 0, asIs = FALSE) {
    
    if (!asIs) {
      obj <- makeBaseCalls(obj)
    }
    
    refSeq <- DNAString(refSeq)
    basecalls <- basecalldf(obj)
    
    refStrandResult <- determineRefStrand(refSeq, basecalls, trim5, trim3)
    refStrand <- refStrandResult[["refStrand"]]
    pairwiseAligned <- refStrandResult[["pa"]]

    #This picks the longest continuous alignment range
    if(pairwiseAligned@pattern@range@width < 10) {
      stop("Seed length not long enough. Must provide at least 10 bases of good 
           matching sequence.")
    }
    
    basecalls <- accountForStartOffset(refSeq, pairwiseAligned, basecalls, trim5)
    basecalls <- updatePrimarySecondary(basecalls)
    
    obj <- updateSangerSeqObj(obj, refStrand, basecalls)
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

updateFinalObj <- function(obj, newObj) {
    if (identical("obj", "objUpdate")) {
      stop("No Change", obj, call. = FALSE)
    }
    else if (exists(obj, envir = "sangerseq", inherits = FALSE)) {
      assign(obj, newObj, envir = "sangerseq")
    } 
    else {
      rebind(obj, newObj, parent.env("sangerseq"))
    }
  }
