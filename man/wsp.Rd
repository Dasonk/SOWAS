\name{wsp}
\alias{wsp}
\title{Wavelet Sample Spectrum}
\usage{
  wsp(ts, s0 = 1, noctave = 5, nvoice = 10, w0 = 2 * pi,
    sw = 0, tw = 0, swabs = 0, siglevel = 0.95,
    critval = NULL, nreal = 1000, arealsiglevel = 0.9,
    kernel = 0, markt = -999, marks = -999,
    logscale = FALSE, plot = TRUE, units = "",
    device = "screen", file = "wsp", color = TRUE,
    pwidth = 10, pheight = 7, labsc = 1, labtext = "",
    sigplot = 3)
}
\arguments{
  \item{ts}{time series object to be transformed}

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

  \item{siglevel}{vector of significance levels for
  pointwise test, e.g. c(0.9,0.95), default 0.95.}

  \item{critval}{matrix of critical values calculated e.g.
  by critical.valuesWSP; each row contains critical values
  for the corresponding chosen significance level.}

  \item{nreal}{number of realizations to estimate critical
  values for the corresponding significance values, default
  1000}

  \item{arealsiglevel}{significance level of the areawise
  test; currently only for siglevel=0.9 and for
  arealsiglevel=0.9 possible, i.e. 90 percent of the area
  of false positive patches is sorted out}

  \item{kernel}{bitmap of the reproducing kernel; if not
  provided, it will be calculated during the areawise test}

  \item{markt}{vector of times to be marked by vertical
  dotted lines; when set to -999 (default), no lines are
  plotted.}

  \item{marks}{vector of scales to be marked by horizontal
  dotted lines; when set to -999 (default), no lines are
  plotted.}

  \item{logscale}{when TRUE, the contours are plotted in
  logarithmic scale}

  \item{plot}{TRUE when graphical output desired}

  \item{units}{character string giving units of the data
  sets. Default: ""}

  \item{device}{"screen" or "ps"}

  \item{file}{character string giving filename of graphical
  output without extension}

  \item{color}{TRUE (default): color plot, FALSE: gray
  scale}

  \item{pwidth}{width of plot in cm}

  \item{pheight}{height of plot in cm}

  \item{labsc}{scale of labels, default: 1, for two-column
  manuscripts: 1.5, for presentations: >2}

  \item{labtext}{puts a label in upper left corner of the
  plot}

  \item{sigplot}{0: no significance test plotted, 1:
  results from pointwise test, 2: results from areawise
  test, 3: results from both tests}
}
\value{
  modulus here: wavelet sample spectrum of dimension
  [length(ts)]x[nvoice*noctave+1] phase not used s0 lowest
  calculated scale in units of the time series noctave
  number of octaves nvoice number of voices per octave w0
  time/frequency resolution omega_0 time vector of times of
  length(intersection of ts1 and ts2) scales vector of
  scales of length nvoice*noctave+1 critval matrix of
  critical values at time/scale matrix of areawise
  significant patches kernel bitmap of reproducing kernel
}
\description{
  This funtion estimates the wavelet spectrum of a time
  series object with the Morlet wavelet.
}
\details{
  A pointwise significance test is performed either by
  providing critical values (e.g. from critical.valuesWSP)
  or by Monte Carlo simulations. If critval and either
  nreal or siglevel is zero, no significance test is
  performed.  If arealsiglevel=0.9, an areawise
  significance test is performed, that sorts out 90 percent
  of the area of false positive patches.  The cone of
  influence is marked by black lines. Values outside the
  cone of influence should be interpreted very carefully,
  as they result from a significant contribution of zero
  padding at the beginning and the end of the time series.
  For a better visualization, additional dotted lines
  marking distinct times or scales might be plotted by
  providing the vectors markt and marks.  In the linear
  plot, values are normalized such that the highest spetral
  power equals one. In the logarithmic plot, values are
  normalized, that the lowest value equals 10^0=1.  The
  function returns an object of type "wt", that might be
  directly plotted by the plot function.
}
\examples{
##
data(nino3)

wspnino3 <- wsp(nino3,s0=0.5,noctave=5,nvoice=10,nreal=100,units="years",arealsiglevel=0)
}
\author{
  D. Maraun
}
\references{
  D. Maraun and J. Kurths, Nonlin. Proc. Geophys. 11:
  505-514, 2004
}
\seealso{
  \code{\link{cwt.ts}}, \code{\link{wcs}},
  \code{\link{wco}}
}
\keyword{power}
\keyword{spectrum}
\keyword{wavelet}

