\name{plotwt}
\alias{plotwt}
\title{Plots wavelet transform}
\usage{
  plotwt(wt, phs, t, s0, noctave, w0, cv = NULL, at = NULL,
    markt = -999, marks = -999, t1 = NULL, t2 = NULL,
    logscale = FALSE, phase = FALSE, units = "",
    plottitle = "Wavelet Plot", device = "screen",
    file = "wavelet", split = FALSE, color = TRUE,
    pwidth = 10, pheight = 7, labsc = 1, labtext = "",
    sigplot = 1, xax = NULL, xlab = NULL, yax = NULL,
    ylab = NULL)
}
\arguments{
  \item{wt}{matrix of real values (e.g. modulus, power,
  coherency) of dimension [length(time
  vector)]x[nvoice*noctave]}

  \item{phs}{matrix of phase values (i.e. [-pi,pi]) of
  dimension [length(time vector)]x[nvoice*noctave]}

  \item{t}{time vector}

  \item{s0}{lowest calculated scale in units of the time
  series}

  \item{noctave}{numbers of octaves}

  \item{w0}{time/frequency resolution omega_0}

  \item{cv}{vector of critical values for each scale of
  length nvoice*noctave+1. If cv=0 is chosen, no critical
  values are plotted.}

  \item{at}{matrix with results from areawise test}

  \item{markt}{vector of times to be marked by vertical
  dotted lines; when set to -999 (default), no lines are
  plotted.}

  \item{marks}{vector of scales to be marked by horizontal
  dotted lines; when set to -999 (default), no lines are
  plotted.}

  \item{t1}{Starting time of the plot. If NULL (default),
  then the plot covers the entire range available}

  \item{t2}{Ending time of the plot. If NULL (default),
  then the plot covers the entire range available}

  \item{logscale}{when TRUE, the contours are plotted in
  logarithmic scale}

  \item{phase}{TRUE if phase is to be plotted}

  \item{units}{character string giving units of the data
  sets; default: ""}

  \item{plottitle}{character string giving plot title}

  \item{device}{"screen" or "ps"}

  \item{file}{character string giving filename of graphical
  output without extension}

  \item{split}{when TRUE, modulus and phase are splitted
  into two files. Default: FALSE}

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
  No value returned
}
\description{
  If plotting results from wsp,wcs or wco, better use
  plot.wt(), which is a wrapper function that calls this
  function.  This function plots results from the functions
  cwt.ts,wsp,wcs,wco.  Modulus and (if possible) phase are
  plotted.
}
\details{
  When cv=0 or -1 is chosen, no critical values are
  plotted. The cone of influence is marked by black lines.
  If cv is a vector, the routine assumes that every value
  represents a critical value, which is constant in scale.
  If cv is a matrix, it assumes every row to contain scale
  dependent critical values, each row stands for one
  significance level.  Values outside the cone of influence
  should be interpreted very carefully, as they result from
  a significant contribution of zero padding at the
  beginning and the end of the time series. For a better
  visualization, additional dotted lines marking distinct
  times or scales might be plotted by providing the vectors
  markt and marks.
}
\examples{
##
}
\author{
  D. Maraun
}
\seealso{
  \code{\link{cwt.ts}}, \code{\link{wsp}},
  \code{\link{wcs}}, \code{\link{wco}}
}
\keyword{plot}
\keyword{transform}
\keyword{wavelet}

