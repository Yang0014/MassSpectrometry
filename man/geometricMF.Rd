\name{geometricMF}
\alias{geometricMF}

\title{
  Euclidean geometric distance matching factor
}
\description{
  The similarity measurement, Euclidean geometric distance matching factor,
  for two spectra.
}
\usage{
  geometricMF(x, y)
}
  
\arguments{
  \item{x}{
    \code{numeric}(n): the peak intensity values of the first spectrum.
  }
  \item{y}{
    \code{numeric}(n): the peak intensity values of the second spectrum,
    with same length of \code{x}.
  }
}
\details{
  This metric measures the geometric distance between two spectra.
  The component of each spectrum is normalised into the unit length.
  Then the normalised vector of spectrum can be considered as a single point 
  a sphere with unit radius in a hyperspace of \code{n} dimensions, 
  where \code{n} is the number of components of the vector.
  Two closer spectra will result in smaller geometric distance.
  
  The inverse of the geometric distance plus 1 is returned in this function,
  as a measure of the similarity of two spectra.
  
  This measurement is one to one correlated with cosine similarity.
}
\value{
  \code{numeric}(1): the similarity score of two spectra between 0.5 and 1.
  1 means perfect identification and 0.5 means the most dissimilarity.
}
\note{
  This similarity measure potentially can be problematic.
  For instance, in the example, a and b are closer to each other than a to c.
  However, considering the magnitude of change of first and last peak, 
  a to b is more dramatically different than a to c.
}
\references{
  Alfassi, Z.B. (2004). On the normalization of a mass spectrum for comparison of two spectra. Journal of the American Society for Mass Spectrometry 15, 385-387.
}
\author{
  Yang Yang
}

\examples{
  a <- c(1, 10, 5, 8)
  b <- c(2, 10, 5, 8)
  c <- c(1, 10, 5, 9)
  geometricMF(a, b)
  geometricMF(a, c)
  geometricMF(b, c)
}