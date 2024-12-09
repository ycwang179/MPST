#' Suppress NOTES for Undefined Global Variables 
#'
#' Certain global variables are used internally in functions such as \code{plot.slice.mpst}
#' for interactive visualization. These variables are dynamically created during execution
#' (e.g., using \code{manipulate::slider}), and therefore static code analysis tools
#' may flag them as undefined. This declaration informs R's static analysis tools
#' (e.g., \code{R CMD check}) that these variables are intentional and should not
#' generate NOTES during package checking.
#'
#' @name globalVariables_MPST
#' @keywords internal
#' @note This declaration suppresses warnings for the following global variables:
#' \itemize{
#'   \item \code{axial_slice}
#'   \item \code{coronal_slice}
#'   \item \code{sagittal_slice}
#'   \item \code{color}
#' }
#' @importFrom utils globalVariables

utils::globalVariables(c("axial_slice", "coronal_slice", "sagittal_slice", "color"))
