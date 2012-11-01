\name{plot.wt}
\alias{plot.wt}
\title{Plots wavelet transform}
\usage{
  plot.wt(wt, markt = -999, marks = -999, t1 = NULL,
    t2 = NULL, logscale = FALSE, phase = FALSE, units = "",
    plottitle = "", device = "screen", file = "wt",
    split = FALSE, color = TRUE, pwidth = 10, pheight = 5,
    labsc = 1.5, labtext = "", sigplot = 3, xax = NULL,
    xlab = NULL, yax = NULL, ylab = NULL)
}
\arguments{
  \item{wt}{Wavelet transformation object}

  \item{markt}{Vector of times to be marked by vertical
  dotted lines. When set to -999 (default), no lines are
  plotted.}

  \item{marks}{Vector of scales to be marked by horizontal
  dotted lines. When set to -999 (default), no lines are
  plotted.}

  \item{t1}{Starting time of the plot. If NULL (default),
  then the plot covers the entire range available}

  \item{t2}{Ending time of the plot. If NULL (default),
  then the plot covers the entire range available}

  \item{logscale}{When TRUE, the contours are plotted in
  logarithmic scale}

  \item{phase}{TRUE if phase is to be plotted}

  \item{units}{character string giving units of the data
  sets. Default: ""}

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
  This function plots objects of type "wt" (wavelet
  transform objects) from the functions wsp,wcs,wco.
  Modulus and (if possible) phase are plotted.
}
\details{
  This function has to be called with an object of type
  "wt" (wavelet transformation object, output from
  wsp,wcs,wco). It calls plotwt().  When cv=0 or -1 is
  chosen, no critical values are plotted. The cone of
  influence is marked by black lines. If cv is a vector,
  the routine assumes that every value represents a
  critical value, which is constant in scale. If cv is a
  matrix, it assumes every row to contain scale dependent
  critical values, each row stands for one significance
  level.  Values outside the cone of influence should be
  interpreted very carefully, as they result from a
  significant contribution of zero padding at the beginning
  and the end of the time series. For a better
  visualization, additional dotted lines marking distinct
  times or scales might be plotted by providing the vectors
  markt and marks.
}
\examples{
##
data(nino3)
nino3wsp <- wsp(nino3,nreal=0)
plot(nino3wsp)
}
\author{
  D. Maraun
}
\seealso{
  \code{\link{plotwt}}, \code{\link{wsp}},
  \code{\link{wcs}}, \code{\link{wco}}
}
\keyword{plot}
\keyword{transform}
\keyword{wavelet}

