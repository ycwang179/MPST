% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.MPST.R
\name{initialize.grid}
\alias{initialize.grid}
\title{Initialize Grid for MPST Plotting}
\usage{
initialize.grid(mfit, Zgrid = NULL, n1 = 101, n2 = 101, n3 = 101)
}
\arguments{
\item{mfit}{An MPST model fit object.}

\item{Zgrid}{An optional grid. If NULL, a grid is created automatically.}

\item{n1}{Number of points in the first dimension.}

\item{n2}{Number of points in the second dimension.}

\item{n3}{Number of points in the third dimension (if applicable).}
}
\value{
A list containing the generated grid and corresponding coordinate vectors.
}
\description{
Generates a grid for plotting based on the input MPST model and dimensionality.
This function is used internally and is not intended for direct use by end users.
}
\keyword{internal}
