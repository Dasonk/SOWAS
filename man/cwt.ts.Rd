\name{cwt.ts}
\alias{cwt.ts}
\title{Continuous Wavelet transformation of time series object}
\usage{
  cwt.ts(ts, s0, noctave = 5, nvoice = 10, w0 = 2 * pi)
}
\arguments{
  \item{ts}{time series object to be transformed}

  \item{s0}{lowest calculated scale in units of the time
  series}

  \item{noctave}{number of octaves}

  \item{nvoice}{number of voices per octave}

  \item{w0}{time/frequency resolution omega_0}
}
\value{
  A matrix containing the (in general complex) Wavelet
  coefficients of dimension [length(intersection of ts1 and
  ts2)]x[nvoice*noctave+1]
}
\description{
  This function calculates the continuous Wavelet
  transformation of a time series object using the Morlet
  Wavelet.
}
\details{
  This function calls the function cwt of the Rwave package
  by Rene Carmona et al. and normalizes the resulting
  voices by sqrt(scale).
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
  \code{\link{wsp}}, \code{\link{wcsp}}, \code{\link{wcoh}}
}
\keyword{continuous}
\keyword{transformation}
\keyword{wavelet}

