getCleanedBasecalls <- function(obj, showtrim, trim5, trim3, averagePosition) {
  
  #Sometimes there are more basecalls than peaks, so the basecalls need to be trimmed
  basecalls1 <- unlist(strsplit(toString(obj@primarySeq), ""))
  basecalls2 <- unlist(strsplit(toString(obj@secondarySeq), ""))
  
  basecalls1 <- basecalls1[1:length(averagePosition)] 
  basecalls2 <- basecalls2[1:length(averagePosition)] 
  
  if(showtrim == FALSE) {
    basecalls1 <- removeTrim(basecalls1, trim5, trim3)
    basecalls2 <- removeTrim(basecalls2, trim5, trim3)
    averagePosition <- averagePosition[(1 + trim5):(length(averagePosition) - trim3)] 
  }
  
  return(c(basecalls1, basecalls2))
}

removeTrim <- function(basecalls, trim5, trim3) {
  
  if(trim5+trim3 > length(basecalls)) {
    basecalls <- ""
  } else {
    basecalls <- basecalls[(1 + trim5):(length(basecalls) - trim3)]
  }
  return(basecalls)
}
