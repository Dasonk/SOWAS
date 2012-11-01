\name{wcs}
\alias{wcs}
\title{Wavelet Sample Cross Spectrum}
\usage{
  wcs(ts1, ts2, s0 = 1, noctave = 5, nvoice = 10,
    w0 = 2 * pi, sw = 0, tw = 0, swabs = 0, markt = -999,
    marks = -999, logscale = FALSE, phase = TRUE,
    plot = TRUE, units = "", device = "screen",
    file = "wcsp", split = FALSE, color = TRUE,
    pwidth = 10, pheight = 7, labsc = 1, labtext = "")
}
\arguments{
  \item{ts1}{first time series object to be transformed}

  \item{ts2}{second time series object to be transformed}

  \item{s0}{lowest calculated scale in units of the time
  series}

  \item{noctave}{number of octaves}

  \item{nvoice}{number of voices per octave}

  \item{w0}{time/frequency resolution omega_0}

  \item{sw}{length of smoothing window in scale direction
  is 2*sw*nvoice+1}

  \item{tw}{length of smoothing window in time direction is
  2*s*tw+1}

  \item{swabs}{length of smoothing window in scale
  direction at scale s is 2*swabs+1}

  \item{markt}{vector of times to be marked by vertical
  dotted lines; when set to -999 (default), no lines are
  plotted.}

  \item{marks}{vector of scales to be marked by horizontal
  dotted lines; when set to -999 (default), no lines are
  plotted.}

  \item{logscale}{when TRUE, the contours are plotted in
  logarithmic scale}

  \item{phase}{TRUE when phase calculation desired}

  \item{plot}{TRUE when graphical output desired}

  \item{units}{character string giving units of the data
  sets. Default: ""}

  \item{device}{"screen" or "ps"}

  \item{file}{character string giving filename of graphical
  output without extension}

  \item{split}{when TRUE, modulus and phase are splitted
  into two files; default: FALSE}

  \item{color}{TRUE (default): color plot, FALSE: gray
  scale}

  \item{pwidth}{width of plot in cm}

  \item{pheight}{height of plot in cm}

  \item{labsc}{scale of labels, default: 1, for two-column
  manuscripts: 1.5, for presentations: >2}

  \item{labtext}{puts a label in upper left corner of the
  plot}
}
\value{
  modulus matrix of modulus of wavelet sample cross
  spectrum of dimension [length(intersection of ts1 and
  ts2)]x[nvoice*noctave+1] phase matrix of phase of wavelet
  sample cross spectrum, same dimension as modulus s0
  lowest calculated scale in units of the time series
  noctave number of octaves nvoice number of voices per
  octave w0 time/frequency resolution omega_0 time vector
  of times of length(intersection of ts1 and ts2) scales
  vector of scales of length nvoice*noctave+1 critval not
  used at not used kernel not used
}
\description{
  This funtion estimates the wavelet cross spectrum of two
  time series objects with the Morlet wavelet.
}
\details{
  \strong{WARNING!} Better do not use this function because
  it is in general easily misinterpreted! A peak in the
  wavelet cross sample spectrum appears in the three cases,
  that either the first processes exhibits a peak, or the
  second process or both. But it does not tell, what case
  is observed. \strong{So in general, a peak in the wavelet
  cross sample spectrum does not imply that the two
  underlying processes are related in any way.} The
  function returns an object of type "wt", that might be
  directly plotted by the plot function.
}
\examples{
##
data(nao)
data(nino3)

# wcs mimics peaks of coherent power, where in reality are non to be
# found, as wco shows (see FAQs on my homepage)
# Thus, never use wcs! :-)
wcsp.nao.nino3 <- wcs(nao,nino3,s0=0.5,noctave=5,nvoice=10)
wcoh.nao.nino3 <- wco(nao,nino3,s0=0.5,noctave=5,nvoice=10,sw=0.5,arealsiglevel=0)
}
\author{
  D. Maraun
}
\references{
  D. Maraun and J. Kurths, Nonlin. Proc. Geophys. 11:
  505-514, 2004
}
\seealso{
  \code{\link{cwt.ts}}, \code{\link{wsp}},
  \code{\link{wco}}
}
\keyword{cross}
\keyword{spectrum}
\keyword{wavelet}

