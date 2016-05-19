### -----------------------------------------------------------------
### normalisation of the mass spectrum
### Exported!
normaliseSpectrum <- function(x, method=c("sum", "max", "unit")){
  method <- match.arg(method)
  if(method == "sum"){
    x <- x / sum(x)
  }else if(method == "max"){
    x <- x / max(x)
  }else{
    x <- x / sqrt(sum(x^2))
  }
  return(x)
}
