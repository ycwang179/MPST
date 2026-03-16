#' Training datasets for the 2D horseshoe example in MPST
#'
#' Simulated training datasets for the 2D horseshoe-domain example used in the
#' simulation studies of \pkg{MPST}. The example is based on the horseshoe
#' domain triangulated into 84 triangles and the target function implemented by
#' \code{func = 15}.
#'
#' The object contains two training datasets corresponding to two noise levels:
#' \describe{
#'   \item{sigma01}{Training data with noise level \eqn{\sigma = 0.1}.}
#'   \item{sigma05}{Training data with noise level \eqn{\sigma = 0.5}.}
#' }
#'
#' Each component is a list with the following entries:
#' \describe{
#'   \item{Y}{A numeric response vector.}
#'   \item{Z}{A numeric matrix of observation locations with 2 columns.}
#'   \item{V}{A numeric matrix of vertex coordinates with 2 columns.}
#'   \item{Tr}{An integer matrix of triangle indices with 3 columns.}
#'   \item{lambda}{The candidate smoothing parameter grid.}
#'   \item{func}{The function identifier used in simulation.}
#'   \item{sigma}{The noise level used in simulation.}
#'   \item{iter}{The random seed used to generate the dataset.}
#'   \item{n}{The sample size.}
#' }
#'
#' @format A named list containing \code{sigma01} and \code{sigma05}.
#' @usage data(ex_hs2d_train)
#' @keywords datasets
"ex_hs2d_train"

#' Prediction datasets for the 2D horseshoe example in MPST
#'
#' Prediction datasets for the 2D horseshoe-domain example used in the
#' simulation studies of \pkg{MPST}. These datasets provide fixed evaluation
#' locations for prediction under the same horseshoe triangulation and target
#' function setting as \code{ex_hs2d_train}.
#'
#' The object contains two prediction datasets corresponding to two noise levels:
#' \describe{
#'   \item{sigma01}{Prediction data with noise level \eqn{\sigma = 0.1}.}
#'   \item{sigma05}{Prediction data with noise level \eqn{\sigma = 0.5}.}
#' }
#'
#' Each component is a list with the following entries:
#' \describe{
#'   \item{Y.grid}{A numeric vector of generated responses on the evaluation set.}
#'   \item{mu.grid}{The underlying mean function values on the evaluation set.}
#'   \item{ind.grid}{An indicator of whether an evaluation point lies inside the domain.}
#'   \item{Z.grid}{A numeric matrix of prediction locations with 2 columns.}
#' }
#'
#' @format A named list containing \code{sigma01} and \code{sigma05}.
#' @usage data(ex_hs2d_pred)
#' @keywords datasets
"ex_hs2d_pred"