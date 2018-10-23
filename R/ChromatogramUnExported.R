cleanBasecalls <- function(basecalls, showtrim, trim5, trim3, averagePosition) {
  
  #Sometimes there are more basecalls than peaks, so the basecalls need to be trimmed
  basecalls <- basecalls[1:length(averagePosition)] 
  if(showtrim == FALSE) {
    basecalls <- removeTrim(basecalls, trim5, trim3)
  }
  return(basecalls)
}

removeTrim <- function(basecalls, trim5, trim3) {
  
  if(trim5+trim3 > length(basecalls)) {
    basecalls <- ""
  } else {
    basecalls <- basecalls[(1 + trim5):(length(basecalls) - trim3)]
  }
  return(basecalls)
}
