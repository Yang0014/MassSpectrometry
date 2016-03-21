
plotPseudoSpectrum = function(x, y, labels=NULL){
  stopifnot(length(x) == length(y))
  #savefig(file=file, width=20, height=10)
  #on.exit(dev.off())
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

plotPseudoGaussianSpectrum = function(x, y, xlim, sd, labels=NULL){
  stopifnot(length(x) == length(y))
  plot(NA,NA,
       xlim=xlim,
       ylim=c(0, 100),
       yaxs="i", bty="l",
       type="n", xlab="mass", ylab="Relative Abundance")
  plotWidths = diff(x)
  plotWidths = c(plotWidths[1], plotWidths, plotWidths[length(plotWidths)])
  for(i in 1:length(x)){
    plot_x = seq(x[i] - plotWidths[i]/2,
            x[i] + plotWidths[i+1]/2, by=0.1)
    plot_y = dnorm(plot_x, mean=x[i], sd=sd)
    plot_y = y[i] / max(plot_y) * plot_y
    lines(plot_x, plot_y)
  }
}

### -----------------------------------------------------------------
### generatePseudoGaussianSpectrum: given x-coordinates of mass or m/z, 
### and the y of peak intensities, it generates the points of gaussian distributsions
### around the peaks for ploting the spectrum.
### Exported!
generatePseudoGaussianSpectrum <- function(x, y, sd=5L, xlim=range(x), step=1L){
  stopifnot(length(x) == length(y))
  plot_x <- seq(xlim[1], xlim[2], by=step)
  ans_y <- numeric(length(plot_x))
  for(i in 1:length(x)){
    plot_y <- dnorm(plot_x, mean=x[i], sd=sd)
    plot_y <- y[i] / max(plot_y) * plot_y
    ans_y <- pmax(ans_y, plot_y, na.rm=TRUE)
  }
  return(list(x=plot_x, y=ans_y))
}



