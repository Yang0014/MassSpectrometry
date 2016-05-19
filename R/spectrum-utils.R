### -----------------------------------------------------------------
### normalisation of the mass spectrum
### Exported!
normaliseSpectrum <- function(x, method=c("sum", "max", "unit")){
  if(any(x < 0)){
    stop("The spectrum intensity values must be non-negative.")
  }
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

### -----------------------------------------------------------------
### medthos of comparison of two spectra
### Exported!

### Euclidean geometric distance matching factor
geometricMF <- function(x, y){
  if(length(x) != length(y)){
    stop("The length of two spectra must be same!")
  }
  x <- normaliseSpectrum(x, method="unit")
  y <- normaliseSpectrum(y, method="unit")
  ans <- 1 + sum((x-y)^2)
  ans <- 1 / ans
  return(ans)
}