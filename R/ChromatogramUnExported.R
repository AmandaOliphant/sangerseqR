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

alignBaseCalls <- function(obj, showcalls, basecallwarning1, basecallwarning2) {
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
  return(obj)
}

