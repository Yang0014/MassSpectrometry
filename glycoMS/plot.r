
plotPseudoSpectrum = function(x, y, labels=NULL, file="output"){
  stopifnot(length(x) == length(y))
  require(monash)
  savefig(file=file, width=20, height=10)
  on.exit(dev.off())
  plot(NA, NA, xlim=range(x), ylim=c(0, 105), type="n", bty="l", 
       yaxs="i", 
       xlab="mass", ylab="Relative abundance")
  for(i in 1:length(x)){
    xToPlot = c(x[i], x[i])
    yToPlot = c(0, y[i])
    lines(xToPlot, yToPlot)
  }
  if(!is.null(labels)){
    text(x=x, y=y+1, labels=labels, cex=0.7)
  }
}


