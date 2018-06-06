removeTrim <- function(basecalls1, basecalls2, trim5, trim3) {
  
  if(trim5+trim3 > length(basecalls1)) {
    basecalls1 <- ""
  }
  else {
    basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
  }
  
  if(trim5+trim3 > length(basecalls2)) {
    basecalls2 <- ""
  }
  else {
    basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
  }
  
  return(list(basecalls1 = basecalls1, basecalls2 = basecalls2))
}

makeSettings <- function(basecalls1, basecalls2, trim5, trim3, 
                         averagePosition, traces, ylim, width) {
  
  indexes <- 1:length(basecalls1)
  trimmed <- indexes <= trim5 | indexes > (length(basecalls1) - trim3) # all false if not trimmed
  
  if (!is.null(trim3)) {
    traces <- traces[1:(min(max(averagePosition, na.rm=TRUE) + 10, 
                            nrow(traces))), ]
  }
  if (!is.null(trim5)) {
    offset <- max(c(1, averagePosition[1] - 10))
    traces <- traces[offset:nrow(traces),]
    averagePosition <- averagePosition - (offset-1)
  }
  
  maxSignal <- apply(traces, 1, max)
  ylims <- c(0, quantile(maxSignal, .75)+ylim*IQR(maxSignal))           
  p <- c(0, averagePosition, nrow(traces))
  midp <- diff(p)/2
  starts <- averagePosition - midp[1:(length(midp)-1)]
  startHets <- starts
  startHets[basecalls1 == basecalls2] <- NA
  ends <- averagePosition + midp[2:(length(midp))]
  endHets <- ends
  endHets[basecalls1 == basecalls2] <- NA
  startTrims <- starts
  startTrims[!trimmed] <- NA
  endTrims <- ends
  endTrims[!trimmed] <- NA
  
  colorTranslate <- c(A="green", C="blue", G="black", T="red")
  colorVector1 <- unname(colorTranslate[basecalls1])
  colorVector1[is.na(colorVector1)] <- "purple"
  colorVector2 <- unname(colorTranslate[basecalls2])
  colorVector2[is.na(colorVector2)] <- "purple"
  
  valuesPerBase <- nrow(traces)/length(basecalls1)
  traceWidth <- width*valuesPerBase
  breaks <- seq(1,nrow(traces), by=traceWidth)
  numPlots <- length(breaks)
  
  return(list(traces = traces, basecalls1 = basecalls1, 
              basecalls2 = basecalls2, averagePosition = averagePosition, ylims = ylims,
              startTrims = startTrims, endTrims = endTrims, startHets = startHets,
              endHets = endHets, colorVector1 = colorVector1, colorVector2 = colorVector2,
              traceWidth = traceWidth, breaks = breaks, numPlots = numPlots))
}

