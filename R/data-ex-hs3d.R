#' Training dataset for the 3D horseshoe example in MPST
#'
#' A simulated training dataset for the 3D irregular horseshoe-domain example
#' used in the visualization studies of \pkg{MPST}. The domain is represented
#' by a tetrahedral partition with 504 tetrahedra.
#'
#' The dataset corresponds to a single representative setting with
#' \code{func = 1}, noise level \eqn{\sigma = 0.1}, and sample size
#' \eqn{n = 10000}.
#'
#' The object is a list with the following entries:
#' \describe{
#'   \item{Y}{A numeric response vector.}
#'   \item{Z}{A numeric matrix of observation locations with 3 columns.}
#'   \item{V}{A numeric matrix of vertex coordinates with 3 columns.}
#'   \item{Tr}{An integer matrix of tetrahedral indices with 4 columns.}
#'   \item{lambda}{The candidate smoothing parameter grid.}
#'   \item{func}{The function identifier used in simulation.}
#'   \item{sigma}{The noise level used in simulation.}
#'   \item{iter}{The random seed used to generate the dataset.}
#'   \item{n}{The sample size.}
#' }
#'
#' @format A list containing training data for the 3D horseshoe example.
#' @usage data(ex_hs3d_train)
#' @keywords datasets
"ex_hs3d_train"

#' Prediction dataset for the 3D horseshoe example in MPST
#'
#' A prediction dataset for the 3D irregular horseshoe-domain example used in
#' the visualization studies of \pkg{MPST}. The domain is represented by a
#' tetrahedral partition with 504 tetrahedra.
#'
#' The object is a list with the following entries:
#' \describe{
#'   \item{Y.grid}{A numeric vector of generated responses on the prediction grid.}
#'   \item{mu.grid}{The underlying mean function values on the prediction grid.}
#'   \item{ind.grid}{An indicator of whether a prediction point lies inside the domain.}
#'   \item{Z.grid}{A numeric matrix of prediction locations with 3 columns.}
#' }
#'
#' @format A list containing prediction data for the 3D horseshoe example.
#' @usage data(ex_hs3d_pred)
#' @keywords datasets
"ex_hs3d_pred"