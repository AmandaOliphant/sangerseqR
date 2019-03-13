trimBasecalls <- function(basecalls, showtrim, trim5, trim3, averagePosition) {
  
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

chromatogramColor <- function(obj, colorTranslate, colorVector1, colorVector2) {
  colorTranslate <- c(A="green", C="blue", G="black", T="red")
  colorVector1 <- unname(colorTranslate[basecalls1])
  colorVector1[is.na(colorVector1)] <- "purple"
  colorVector2 <- unname(colorTranslate[basecalls2])
  colorVector2[is.na(colorVector2)] <- "purple"
  return(obj)
}
