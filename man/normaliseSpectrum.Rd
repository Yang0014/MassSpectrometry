\name{normaliseSpectrum}
\alias{normaliseSpectrum}

\title{
  Spectrum normalisation
}
\description{
  Three normalisation methods for spectrum.
}
\usage{
  normaliseSpectrum(x, method = c("sum", "max", "unit"))
}
\arguments{
  \item{x}{
    \code{numeric}(n): the intensity values at each \emph{m/z}.
  }
  \item{method}{
    \code{character}(1): the method for normalisation.
    "sum": the sum of all the intensity values is taken as unity, and each 
    peak is normalised by dividing by the sum.
    "max": each peak is normalised by dividing by the highest peak.
    "unit": each peak is normalised by dividing by the squared root of the sum
    of the squred peaks,
    \eqn{\sqrt{\sum_{j=1}^{n}I_j^2}}{sqrt(sum(I_i^2))}, 
    which leads to each spectrum being a vector of unit length.
  }
}
\details{
  Rasmussen and Isenhour studied these three normalisation methods and
  they found that all the normalization methods and search methods 
  gave similar results, 
  although the "sum" method seems the to be best normalization method.
}
\value{
  \code{numeric}(n) vector of normalised peaks.
}
\references{
  Rasmussen, G.T., and Isenhour, T.L. (1979). The Evaluation of Mass Spectral Search Algorithms. J. Chem. Inf. Comput. Sci. 19, 179-186.
}
\author{
  Ge Tan
}

\examples{
  x <- c(50,100,10,200)
  normaliseSpectrum(x, method="sum")
  normaliseSpectrum(x, method="max")
  normaliseSpectrum(x, method="unit")
}