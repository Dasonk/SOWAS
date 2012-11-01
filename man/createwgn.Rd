\name{createwgn}
\alias{createwgn}
\title{White noise time series}
\usage{
  createwgn(N, sig, dt)
}
\arguments{
  \item{N}{Number of data points}

  \item{sig}{sqrt of variance of the process}

  \item{dt}{sampling time}
}
\value{
  Returns time series object of white noise
}
\description{
  Creates white noise time series object
}
\examples{
##
}
\author{
  D. Maraun
}
\seealso{
  \code{\link[stats]{rnorm}}, \code{\link{create.ar}}
}
\keyword{noise}
\keyword{white}

