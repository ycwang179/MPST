% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPST.R, R/fit.MPST.R
\name{fit}
\alias{fit}
\alias{fit.MPST}
\title{Generic function for model fitting}
\usage{
fit(formula, lambda = NULL, method = NULL, P.func = NULL, data = list())

\method{fit}{MPST}(formula, lambda = NULL, method = NULL, P.func = NULL, data = list())
}
\arguments{
\item{formula}{A formula specifying the model, e.g., `Y ~ m(Z, V, Tr, d, r)`. 
- `Y`: The response variable observed over the domain.
- `Z`: Matrix of observation coordinates (\code{n} by \code{k}). Rows represent points in 
  2D or 3D space (\code{k = 2} or \code{k = 3}). \( k \) is the dimension of the observed 
  points, where \( k = 2 \) for 2D and \( k = 3 \) for 3D.
- `V`: Matrix of vertices (\code{nV} by \code{k}). Rows represent coordinates of vertices 
  in the triangulation.
- `Tr`: Triangulation matrix (\code{nT} by \code{k+1}). Rows represent vertex indices:
  - For 2D: Rows have three indices for triangles.
  - For 3D: Rows have four indices for tetrahedra.
- `d`: Degree of piecewise polynomials (default: \code{5}). \code{-1} represents piecewise constants.
- `r`: Smoothness parameter (default: \code{1}, where \code{0 <= r < d}).}

\item{lambda}{A numeric vector of tuning parameters for regularization. Defaults to 
\eqn{10^(-6,-5.5,-5,...,5,5.5,6)}.}

\item{method}{A character string specifying the learning method. If not specified, defaults to `"G"` (Global learning).
- `"G"`: Global learning.
- `"D"`: Distributed learning.}

\item{P.func}{An integer specifying the parallelization method for distributed learning. Defaults to \code{2}:
- `1`: Use `mclapply`.
- `2`: Use `parLapply`.}

\item{data}{(Optional) A list containing the following components:
- `Y`: The response variable observed over the domain.
- `Z`: Matrix of observation coordinates.
- `V`: Matrix of triangulation vertices.
- `Tr`: Triangulation matrix.
- `d`: Degree of piecewise polynomials.
- `r`: Smoothness parameter.}
}
\value{
An object of class `"MPST"` with the following components:
- `gamma.hat`: Estimated spline coefficients from the fitted model.
- `Y.hat`: Predicted values based on the fitted model.
- `mse`: Mean squared error of the model (computed during fitting).
- `mise`: Mean integrated squared error.
- `method`: The learning method used ("G" or "D").
- `formula`: The formula provided during fitting.
}
\description{
The `fit` function is a generic function for model fitting. It takes a formula, 
penalty parameters, and a method as inputs and delegates the computation 
to specific methods like `fit.MPST()`.

`fit.MPST()` fits a Multivariate Penalized Spline over Triangulation (MPST)
model using global (`"G"`) or distributed (`"D"`) learning methods. The function can extract
required parameters (`Y`, `Z`, `V`, `Tr`, `d`, `r`) either from the formula or from the provided
data list, prioritizing values in the formula if both sources are available.
}
\details{
- This function extracts required components (`Y`, `Z`, `V`, `Tr`, `d`, `r`) from the `formula` 
  using `interpret.mpst()`. If a component is not available in the formula, it falls back to the 
  `data` argument.
- If a required component is missing in both `formula` and `data`, the function raises an error.
- The `method` parameter specifies the learning mode, and `lambda` allows for fine-tuning the 
  regularization.
}
