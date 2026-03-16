#' Training datasets for 3D cube examples in MPST
#'
#' Simulated training datasets for a 3D cube domain under different
#' tetrahedral resolutions. These datasets are included for illustrating
#' model fitting with \code{fit.MPST()}.
#'
#' The object contains three training datasets:
#' \describe{
#'   \item{tet48}{Training data based on a tetrahedral partition with 48 tetrahedra.}
#'   \item{tet72}{Training data based on a tetrahedral partition with 72 tetrahedra.}
#'   \item{tet108}{Training data based on a tetrahedral partition with 108 tetrahedra.}
#' }
#'
#' Each component is a list with the following entries:
#' \describe{
#'   \item{Y}{A numeric response vector.}
#'   \item{Z}{A numeric matrix of observation locations with 3 columns.}
#'   \item{V}{A numeric matrix of vertex coordinates with 3 columns.}
#'   \item{Tr}{An integer matrix of tetrahedral indices with 4 columns.}
#'   \item{mu}{The underlying mean function values used in simulation.}
#' }
#'
#' @format A named list containing \code{tet48}, \code{tet72}, and \code{tet108}.
#' @usage data(ex_cube_train)
#' @keywords datasets
"ex_cube_train"

#' Prediction datasets for 3D cube examples in MPST
#'
#' Prediction grids for 3D cube examples under different tetrahedral
#' resolutions. These datasets are included for illustrating prediction with
#' \code{predict.MPST()}.
#'
#' The object contains three prediction datasets:
#' \describe{
#'   \item{tet48}{Prediction data corresponding to the 48-tetrahedra training example.}
#'   \item{tet72}{Prediction data corresponding to the 72-tetrahedra training example.}
#'   \item{tet108}{Prediction data corresponding to the 108-tetrahedra training example.}
#' }
#'
#' Each component is a list with the following entries:
#' \describe{
#'   \item{Z.grid}{A numeric matrix of prediction locations with 3 columns.}
#'   \item{Y.grid}{A numeric vector of generated responses on the prediction grid.}
#'   \item{mu.grid}{The underlying mean function values on the prediction grid.}
#'   \item{ind.grid}{An indicator of whether a grid point lies inside the domain.}
#' }
#'
#' @format A named list containing \code{tet48}, \code{tet72}, and \code{tet108}.
#' @usage data(ex_cube_pred)
#' @keywords datasets
"ex_cube_pred"