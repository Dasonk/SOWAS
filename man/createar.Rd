\name{createar}
\alias{createar}
\title{AR process time series}
\usage{
  createar(N, a, sig, dt)
}
\arguments{
  \item{N}{Number of data points}

  \item{a}{vector of AR coefficients}

  \item{sig}{sqrt of variance of the process}

  \item{dt}{sampling time}
}
\value{
  Returns time series object of AR process
}
\description{
  Creates AR process time series object
}
\examples{
##
}
\author{
  D. Maraun
}
\seealso{
  \code{\link[stats]{arima.sim}}, \code{\link{create.wgn}}
}
\keyword{AR}
\keyword{process}

