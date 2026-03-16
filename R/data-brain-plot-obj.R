#' Precomputed Brain Slice Volumes for MPST Visualization
#'
#' A collection of precomputed 3D brain volumes for demonstrating
#' slice visualization with \code{\link{plot.MPST}}.
#'
#' The object \code{ex_brain_plot_obj} is a named list. Each component is an
#' MPST-compatible precomputed plotting object that can be directly passed to
#' \code{plot.MPST(..., mview = "slice")}.
#'
#' Available components include:
#' \itemize{
#'   \item \code{global}: precomputed global learning MPST estimate
#'   \item \code{distributed}: precomputed distributed learning MPST estimate
#'   \item \code{tps1}: precomputed thin-plate spline estimate
#'   \item \code{tps2}: precomputed tensor product spline estimate
#'   \item \code{truth}: precomputed true function volume
#'   \item \code{noise}: precomputed noisy volume
#' }
#'
#' Each component is a list with:
#' \itemize{
#'   \item \code{Yarray}: a 3D array of voxel values
#'   \item \code{xgrid}: grid values along the first dimension
#'   \item \code{ygrid}: grid values along the second dimension
#'   \item \code{zgrid}: grid values along the third dimension
#'   \item \code{slice_style}: plotting style, set to \code{"brain"}
#' }
#'
#' @format A named list of six MPST-compatible precomputed plotting objects.
#'
#' @examples
#' data(ex_brain_plot_obj)
#'
#' \dontrun{
#' plot.MPST(ex_brain_plot_obj$global, mview = "slice")
#' plot.MPST(ex_brain_plot_obj$distributed, mview = "slice", slice_style = "brain")
#' plot.MPST(ex_brain_plot_obj$truth, mview = "slice", slice_style = "brain")
#' }
#'
#' @name ex_brain_plot_obj
NULL