\name{generatePseudoGaussianSpectrum}
\alias{generatePseudoGaussianSpectrum}

\title{
  Generate the pseudo Gaussian distribution of peaks
}
\description{
  Given x-coordinates of mass or m/z, and the y of peak intensities, 
  it generates the points of Gaussian distributsions around the peaks for 
  ploting the spectrum.
}
\usage{
  generatePseudoGaussianSpectrum(x, y, sd=5L, xlim=range(x), step = 1L)
}

\arguments{
  \item{x}{
    \code{numeric}(n): the x-coordinates of mass or m/z.
    }
  \item{y}{
    \code{numeric}(n): the y values of peak intensities.
  }
  \item{sd}{
    \code{numeric}(1): the desired standard deviation of the
    Gaussian distribution.
  }
  \item{xlim}{
    \code{numeric}(2): the range for generating points.
    By default, it's the range of x
  }
  \item{step}{
    \code{numeric}(1): the step size of generating points. Smaller step size
    gives smoother lines.
  }
}
\details{
  For each pair of x and y, we generate the Gaussian distribution with \code{x}
  as mean, \code{y} as the height, \code{sd} as the standard deviation.
}
\value{
  A \code{list} of x values and y values for plotting the spectrum.
}

\author{
  Yang Yang
}

\examples{
  x <- c(100, 500, 800)
  y <- c(50, 100, 20)
  ans <- generatePseudoGaussianSpectrum(x, y, sd=5L, xlim=c(1, 1000), step = 1L)
  plot(ans[["x"]], ans[["y"]], type="l", xlab="m/z", ylab="Intensity")
}
