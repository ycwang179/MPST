% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataGenerator3D.R
\name{dataGenerator3D}
\alias{dataGenerator3D}
\title{Data Generator for Numerical Example}
\usage{
dataGenerator3D(Z, V, Tr, func = 1, sigma = 1, seed = 2023)
}
\arguments{
\item{Z}{The coordinates of dimension \code{n} by three. Each row is the 3D coordinates of a point.
\cr}

\item{V}{The \code{N} by three matrix of vertices of a tetrahedron, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
\cr}

\item{Tr}{The tetrahedral partition matrix of dimension \code{nT} by four, where \code{nT} is the number of tetrahedrons in the partition. Each row is the indices of vertices in \code{V}.
\cr}
}
\value{
A data frame.
}
\description{
Data Generator for Numerical Example
}
