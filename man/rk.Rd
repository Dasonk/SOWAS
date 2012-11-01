\name{rk}
\alias{rk}
\title{Reproducing Kernel of the Morlet Wavelet}
\usage{
  rk(N = 1000, s = 8, noctave = 5, nvoice = 10,
    w0 = 2 * pi, plot = TRUE)
}
\arguments{
  \item{N}{length of time series (dt set to one!)}

  \item{s}{scale at which the r.k. is to be calculated}

  \item{noctave}{number of octaves}

  \item{nvoice}{number of voices per octave}

  \item{w0}{time/frequency resolution omega_0}

  \item{plot}{TRUE when grapical output desired}
}
\value{
  Matrix of r.k. of dimension [N]x[nvoice*noctave+1]
}
\description{
  This funtion calculates the reproducing Kernel of the
  Morlet wavelet.
}
\details{
  This function calculates the reproducing kernel of the
  Morlet wavelet at a given scale, i.e. the internal
  correlations of the wavelet coefficients for white noise
  at this scale. These are responsible for the spurious
  patches in wavelet (cross) spectral analysis and are an
  intrinsic property of any time/frequency spectral
  analysis.
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
  \code{\link{cwt.ts}}
}
\keyword{kernel}
\keyword{reproducing}

