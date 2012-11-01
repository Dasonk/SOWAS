\name{criticalvaluesWSP}
\alias{criticalvaluesWSP}
\title{Estimates critical values for Wavelet spectra}
\usage{
  criticalvaluesWSP(ts, s0 = 1, noctave = 5, nvoice = 10,
    w0 = 2 * pi, swabs = 0, tw = 0, siglevel = 0.95,
    nreal = 1000)
}
\arguments{
  \item{ts}{time series object}

  \item{s0}{lowest calculated scale in units of the time
  series}

  \item{noctave}{number of octaves}

  \item{nvoice}{number of voices per octave}

  \item{w0}{time/frequency resolution \omega_0}

  \item{swabs}{length of smoothing window is 2swabs+1}

  \item{tw}{length of smoothing window in time direction is
  2*s*tw+1}

  \item{siglevel}{significance level, e.g. 0.9, 0.95 or
  0.99.  siglevel might also be a vector, e.g. c(0.9,0.95)
  to plot more contourlines.}

  \item{nreal}{number of realizations to estimate critical
  values for the corresponding significance values, default
  1000}
}
\value{
  Returns a matrix of scale DEPENDENT critical values. The
  number of rows of the matrix corresponds to the number of
  chosen significance values. The number of columns equals
  the number of scales.
}
\description{
  This function estimates critical values for the Wavelet
  spectra by means of Monte Carlo simulations. Null
  Hypothesis is an AR1 (red noise) process fitted to the
  given time series.
}
\details{
  nreal might be chosen as 100 for a rough estimate of
  significance. However, it is for sure not suitable to
  reliably distinguish between 95 and 99 percent
  significance values. In this case, at least nreal=1000
  should be chosen.
}
\examples{
##
}
\author{
  D. Maraun
}
\references{
  D. Maraun and J. Kurths, Nonlin. Proc. Geophys. 11:
  505-514, 2004
}
\seealso{
  \code{\link{wcoh}}
}
\keyword{coherency}
\keyword{critival}
\keyword{significance}
\keyword{values}
\keyword{wavelet}

