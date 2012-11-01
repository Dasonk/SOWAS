\name{criticalvaluesWCO}
\alias{criticalvaluesWCO}
\title{Estimates critical values for wavelet coherency}
\usage{
  criticalvaluesWCO(s0, noctave, nvoice, w0, swabs, tw,
    siglevel = 0.95)
}
\arguments{
  \item{s0}{lowest calculated scale in units of the time
  series}

  \item{noctave}{number of octaves}

  \item{nvoice}{number of voices per octave}

  \item{w0}{time/frequency resolution omega_0}

  \item{swabs}{length of smoothing window in scale
  direction at scale s is 2*swabs+1}

  \item{tw}{length of smoothing window in time direction is
  2*s*tw+1}

  \item{siglevel}{significance level, e.g. 0.9, 0.95 or
  0.99.  At the moment, only these values are possible.
  siglevel might also be a vector, e.g. c(0.9,0.95) to plot
  more contourlines.}
}
\value{
  Returns a scale independent critical value. If siglevel
  is a vector of multiple significance values, also the
  return value is a vector of the same length.
}
\description{
  This function estimates critical values for the wavelet
  coherency between two processes using the Morlet wavelet
  based on a formula derived from Monte Carlo simulations.
}
\details{
  The process dependency appeared to be rather marginal.
  Thus we performed MC simulations with two Gaussian White
  Noise processes for the listed significance levels for
  different smoothing windows and time/frequency
  resolutions.
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
\keyword{critical}
\keyword{significance}
\keyword{values}
\keyword{wavelet}

