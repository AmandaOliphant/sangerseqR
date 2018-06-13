removeTrim <- function(basecalls, trim5, trim3) {
  
  if(trim5+trim3 > length(basecalls)) {
    basecalls <- ""
  } else {
    basecalls <- basecalls[(1 + trim5):(length(basecalls) - trim3)]
  }
  return(basecalls)
}
