removeTrim <- function(basecalls1, basecalls2, trim5, trim3) {
  
  if(trim5+trim3 > length(basecalls1)) {
    basecalls1 <- ""
  } else {
    basecalls1 <- basecalls1[(1 + trim5):(length(basecalls1) - trim3)]
  }
  
  if(trim5+trim3 > length(basecalls2)) {
    basecalls2 <- ""
  } else {
    basecalls2 <- basecalls2[(1 + trim5):(length(basecalls2) - trim3)]
  }
  
  return(list(basecalls1 = basecalls1, basecalls2 = basecalls2))
}
