#' Training datasets for 2D square examples in MPST
#'
#' Simulated training datasets for a square domain under different
#' triangulation schemes, surface functions, and noise levels.
#' These datasets are included for illustrating model fitting with
#' \code{fit.MPST()}.
#'
#' The object contains two surface settings:
#' \describe{
#'   \item{m1}{A moderately oscillating surface defined by
#'   \eqn{m_1(z_1,z_2) = -\sin\{3\pi(z_1+0.25)\} + \sin(3\pi z_2)}.}
#'   \item{m2}{A highly oscillating surface defined by
#'   \eqn{m_2(z_1,z_2) = -\sin\{10\pi(z_1+0.25)\} + \sin(10\pi z_2)}.}
#' }
#'
#' For each surface, one noise level is included:
#' \describe{
#'   \item{sigma01}{Datasets generated with \code{sigma = 0.1}.}
#' }
#'
#' For each combination of surface and noise level, six triangulation schemes
#' are provided:
#' \describe{
#'   \item{tri8}{Triangulation with 8 triangles.}
#'   \item{tri18}{Triangulation with 18 triangles.}
#'   \item{tri32}{Triangulation with 32 triangles.}
#'   \item{tri50}{Triangulation with 50 triangles.}
#'   \item{tri72}{Triangulation with 72 triangles.}
#'   \item{tri98}{Triangulation with 98 triangles.}
#' }
#'
#' Each training dataset is a list with the following entries:
#' \describe{
#'   \item{Y}{A numeric response vector.}
#'   \item{Z}{A numeric matrix of observation locations with 2 columns.}
#'   \item{V}{A numeric matrix of vertex coordinates with 2 columns.}
#'   \item{Tr}{An integer matrix of triangle indices with 3 columns.}
#' }
#'
#' These datasets were generated with sample size \code{n = 20000}
#' and random seed \code{2026}.
#'
#' @format A nested list of training datasets indexed by
#' \code{m1}/\code{m2}, \code{sigma01}, and
#' \code{tri8}/\code{tri18}/\code{tri32}/\code{tri50}/\code{tri72}/\code{tri98}.
#' @usage data(ex_square_train)
#' @keywords datasets
"ex_square_train"


#' Prediction datasets for 2D square examples in MPST
#'
#' Prediction grids for 2D square examples under different triangulation
#' schemes, surface functions, and noise levels.
#' These datasets are included for illustrating prediction with
#' \code{predict.MPST()}.
#'
#' The object contains two surface settings:
#' \describe{
#'   \item{m1}{Prediction datasets for the surface
#'   \eqn{m_1(z_1,z_2) = -\sin\{3\pi(z_1+0.25)\} + \sin(3\pi z_2)}.}
#'   \item{m2}{Prediction datasets for the surface
#'   \eqn{m_2(z_1,z_2) = -\sin\{10\pi(z_1+0.25)\} + \sin(10\pi z_2)}.}
#' }
#'
#' For each surface, one noise level is included:
#' \describe{
#'   \item{sigma01}{Prediction datasets generated with \code{sigma = 0.1}.}
#' }
#'
#' For each combination of surface and noise level, six triangulation schemes
#' are provided:
#' \describe{
#'   \item{tri8}{Prediction data corresponding to the 8-triangle training example.}
#'   \item{tri18}{Prediction data corresponding to the 18-triangle training example.}
#'   \item{tri32}{Prediction data corresponding to the 32-triangle training example.}
#'   \item{tri50}{Prediction data corresponding to the 50-triangle training example.}
#'   \item{tri72}{Prediction data corresponding to the 72-triangle training example.}
#'   \item{tri98}{Prediction data corresponding to the 98-triangle training example.}
#' }
#'
#' Each prediction dataset is a list with the following entries:
#' \describe{
#'   \item{Z.grid}{A numeric matrix of prediction locations with 2 columns.}
#'   \item{Y.grid}{A numeric vector of generated responses on the prediction grid.}
#'   \item{mu.grid}{The underlying mean function values on the prediction grid.}
#'   \item{ind.grid}{An indicator for whether a grid point lies inside the domain.}
#' }
#'
#' These prediction datasets were generated with random seed \code{2026}.
#'
#' @format A nested list of prediction datasets indexed by
#' \code{m1}/\code{m2}, \code{sigma01}, and
#' \code{tri8}/\code{tri18}/\code{tri32}/\code{tri50}/\code{tri72}/\code{tri98}.
#' @usage data(ex_square_pred)
#' @keywords datasets
"ex_square_pred"