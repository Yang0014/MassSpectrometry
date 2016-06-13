\name{readMZ}
\alias{readMZ}
\title{
  Read mzXML files
}
\description{
  Read mzXML files.
}
\usage{
readMZ(fns, scanCounts = NULL, starts = NULL, ends = NULL)
}

\arguments{
  \item{fns}{
    \code{character}(n): the filenames of mzXML files.
  }
  \item{scanCounts}{
    \code{integer}(n): the scanCount to use in each file. If NULL,
    the last scan will be used.
  }
  \item{starts}{
    \code{numeric}(n): the start coordinates of the ranges.
  }
  \item{ends}{
    \code{numeric}(n): the end coordinates of the ranges.
  }
}
\details{
  The peak intensities are normalised to the maximal intensity of 100 
  within eahc file.
}
\value{
  When \code{starts} and \code{ends} are NULL, a list of \code{matrix} with
  intensities will be returned.
  When they are given, a \code{data.frame} of the added intensities
  within the ranges will be returned.
}
\author{
  Yang Yang
}

\examples{
  mzFns <- system.file(c("threonine/threonine_i2_e35_pH_tree.mzXML",
                         "lockmass/LockMass_test.mzXML"),
                       package = "msdata")

  ## Read all the mzXML data
  allIntensities <- readMZ(mzFns)

  ## Read the peaks data within certain ragnges
  starts <- c(50, 70, 80)
  ends <- c(55, 75, 85)
  rangedIntensities <- readMZ(mzFns, starts=starts, ends=ends)
}