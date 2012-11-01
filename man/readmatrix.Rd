\name{readmatrix}
\alias{readmatrix}
\title{Loads matrix from file}
\usage{
  readmatrix(file, M)
}
\arguments{
  \item{file}{a character string giving the name of the
  file to load.}

  \item{M}{number of columns of the matrix to be loaded.}
}
\value{
  Returns matrix with M columns
}
\description{
  This function loads a tab separated matrix with M columns
  in ASCII-format into an R matrix.
}
\examples{
##
}
\author{
  D. Maraun
}
\seealso{
  \code{\link{readts}}, \code{\link{writematrix}}
}
\keyword{file}

