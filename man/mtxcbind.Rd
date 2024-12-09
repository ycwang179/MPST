\name{mtxcbind}
\alias{mtxcbind}
\title{Combine Matrices Column-wise}
\usage{
mtxcbind(a, b)
}
\arguments{
  \item{a}{A numeric matrix.}

  \item{b}{A numeric matrix.}
}
\value{
A numeric matrix combining columns of \code{a} and \code{b}.
}
\description{
This function combines two numeric matrices column-wise.
}
\examples{
a <- matrix(1:4, nrow = 2)
b <- matrix(5:8, nrow = 2)
mtxcbind(a, b)
}
\seealso{
Other matrix manipulation functions.
}
