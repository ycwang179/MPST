# MPST.R - Core Functions and Generic Definitions

#' Generic function for model fitting
#'
#' The `fit` function is a generic function for model fitting. It takes a formula, 
#' penalty parameters, and a method as inputs and delegates the computation 
#' to specific methods like `fit.MPST()`.
#'
#' @rdname fit
#' @export
fit <- function(formula, lambda = NULL, method = NULL, P.func = NULL, data = list()) {
  class(formula) <- c(class(formula), "MPST")
  if ("MPST" %in% class(formula) || inherits(data, "MPST")) {
    return(fit.MPST(formula, lambda, method, P.func, data))
  } else {
    # If no MPST compatibility, fallback to generic method
    UseMethod("fit")
  }
}

#' Generic function for model prediction
#'
#' The `predict` function is a generic function for generating predictions from a fitted model object. 
#' Specific methods like `predict.MPST()` provide implementations for specific model types, such as 
#' Multivariate Penalized Spline over Triangulation (MPST).
#'
#' @rdname predict
#' @export
predict <- function(formula, lambda = NULL, method = NULL, P.func = NULL, data = list(), data.pred = list()) {
  class(formula) <- c(class(formula), "MPST")
  # Check if formula and data are MPST-compatible
  if ("MPST" %in% class(formula) || inherits(data, "MPST")) {
    # Check if data.pred is MPST-compatible
    if (inherits(data.pred, "MPST")) {
      # Call predict.MPST if all conditions are met
      return(predict.MPST(formula, lambda, method, P.func, data, data.pred))
    } else {
      stop("'data.pred' must inherit from class 'MPST' for predict.MPST().")
    }
  }
  # If no MPST compatibility, fallback to generic method
  UseMethod("predict")
}

#' Print Method for MPST Models
#'
#' Prints a summary of an MPST model object to the console, including key attributes and results, based on the function type (\code{"fit"}, \code{"predict"}, or \code{"basis"}).
#'
#' @param x An object of class \code{"MPST"}, representing a fitted model, prediction results, or basis functions.
#'   The object must include the following fields:
#'   \itemize{
#'     \item \code{formula}: The formula used in the model.
#'     \item \code{func}: The function type, one of \code{"fit"}, \code{"predict"}, or \code{"basis"}.
#'     \item \code{method}: The method type, either \code{"G"} (Global) or \code{"D"} (Distributed Learning).
#'     \item Additional fields depending on \code{func}, such as \code{gamma.hat}, \code{Y.hat}, \code{mise}, or \code{basis.all}.
#'   }
#' @param ... Additional arguments (currently unused).
#' 
#' @return Invisibly returns the input \code{x} after printing its summary to the console.
#' 
#' @examples
#' # Example MPST object
#' x <- list(
#'   formula = y ~ m(Z, V, Tr, d, r),
#'   func = "fit",
#'   method = "G",
#'   mse = 0.025,
#'   gamma.hat = runif(10),
#'   Y.hat = runif(10),
#'   d = 2
#' )
#' class(x) <- "MPST"
#' print.MPST(x)
#' 
#' @export
print.MPST <- function(x, ...) {
  if (!inherits(x, "MPST")) {
    stop("Object must be of class 'MPST'")
  }
  options(max.print = 10)
  cat("\nFormula: ")
  print(x$formula)
  
  if (x$func == "fit") {
    # handle fit Situations
    if (x$method == "G") {
      cat("\nMethod: Global\n")
    } else if (x$method == "D") {
      cat("\nMethod: Distributed Learning\n")
      if (x$P.func == 1) {
        cat("\nParallelization function: mclapply() with", x$N.cores, "cores\n")
      } else if (x$P.func == 2) {
        cat("\nParallelization function: parLapply() with", x$N.cores, "cores\n")
      }
    }
    cat("\nd:", x$d, "\n")
    cat("\nMSE:", round(x$mse, 4), "\n")  
    cat("\ngamma_hat:\n")
    print(as.vector(x$gamma.hat))
    cat("\ny_hat:\n")
    print(x$Y.hat)
    
  } else if (x$func == "predict") {
    # handle predict Situations
    if (x$method == "G") {
      cat("\nMethod: Global\n")
    } else if (x$method == "D") {
      cat("\nMethod: Distributed Learning\n")
      if (x$P.func == 1) {
        cat("\nParallel function: mclapply() with", x$N.cores, "cores\n")
      } else if (x$P.func == 2) {
        cat("\nParallel function: parLapply() with", x$N.cores, "cores\n")
      }
    }
    cat("\nd:", x$d, "\n")
    cat("\nMISE:", round(x$mise, 4), "\n")  
    cat("\ny_predict:\n")
    print(x$Ypred)
    
  } else if (x$func == "basis") {
    # handle basis Situations
    if (x$method == "G") {
      cat("\nMethod: Global\n")
      cat("\nbasis:\n")
      print(x$B)
    } else if (x$method == "D") {
      cat("\nMethod: Distributed Learning\n")
      for (item in 1:length(x$basis.all)) {
        cat("\nbasis for triangle", item, ":\n")
        print(x$basis.all[[item]]$B.star)
      }
    }
  }
  
  invisible(x)
}

# Define the interpret.mpst function
#' Interpret an MPST Formula
#'
#' @description Parses and interprets an MPST formula for use in model fitting. 
#' This function extracts the response variable, if present, and parses model parameters 
#' provided in the formula.
#'
#' @param mpstf An MPST formula object, such as \code{Y ~ m(Z, V, Tr, d, r)}.
#' - `Y`: The response variable observed over the domain.
#' - `Z`: Matrix of observation coordinates.
#' - `V`: Matrix of triangulation vertices.
#' - `Tr`: Triangulation matrix.
#' - `d`: Degree of piecewise polynomials.
#' - `r`: Smoothness parameter.
#' @param extra.special (Optional) Additional special terms to recognize in the formula. Defaults to \code{NULL}.
#' @return A list containing:
#' - `Y`: The response variable, extracted from the formula or set to \code{NA} if not found.
#' - `Z`: The matrix of observation coordinates.
#' - `V`: The matrix of triangulation vertices.
#' - `Tr`: The triangulation matrix.
#' - `d`: Degree of piecewise polynomials.
#' - `r`: Smoothness parameter.
#' @details
#' - The function first checks for the presence of a response variable in the formula.
#' - If the response variable exists, it is extracted and evaluated in the parent environment.
#' - If the response variable is not found, it is set to \code{NA}.
#' - Model parameters specified within the \code{m()} function are parsed and returned in a list.
#' - The function expects the formula to follow the format \code{Y ~ m(Z, V, Tr, d, r)}.
#'
#' @keywords internal
interpret.mpst <- function(mpstf, extra.special = NULL) {
  # Parse the formula and check for specific terms
  mpst.tf <- terms.formula(mpstf, specials = c("m", extra.special)) 
  mpst.terms <- attr(mpst.tf, "term.labels") 
  
  # Check for the presence of a response variable
  if (attr(mpst.tf, "response") > 0) {  
    response_name <- as.character(attr(mpst.tf, "variables")[2])
    # Attempt to locate the response variable in the parent environment
    if (exists(response_name, envir = parent.frame())) {
      mpst.response <- eval(parse(text = response_name), envir = parent.frame())
    } else {
      mpst.response <- NA 
    }
  } else { 
    mpst.response <- NA
  }
  
  # Parse the parameter section
  mpst.parameters <- eval(parse(text = paste(mpst.terms[1])), envir = parent.frame())
  mpst.parameters$Y <- mpst.response
  
  return(mpst.parameters)
}

# Define the m function
#' Generate MPST Parameters
#'
#' @description Internal function to create a list of parameters for MPST. This function validates the input data and parameters to ensure consistency for either 2D or 3D triangulations.
#' @param ... Additional arguments (not used).
#' @param Y The response variable observed over the domain.
#' @param Z Matrix of observation coordinates (\code{n} by \code{k}). Rows represent points in 2D or 3D space (\code{k = 2} or \code{k = 3}). \( k \) is the dimension of the observed points, where \( k = 2 \) for 2D and \( k = 3 \) for 3D.
#' @param V Matrix of vertices (\code{nV} by \code{k}). Rows represent coordinates of vertices in the triangulation.
#' @param Tr Triangulation matrix (\code{nT} by \code{k+1}). Rows represent vertex indices:
#'   - For 2D: Rows have three indices for triangles.
#'   - For 3D: Rows have four indices for tetrahedra.
#' @param d Degree of piecewise polynomials. Must be a numeric integer value greater than or equal to 1 (default: \code{5}). \code{-1} represents piecewise constants.
#' @param r Smoothness parameter (default: \code{1}). Must satisfy the condition \code{0 <= r < d}.
#' @return A list containing the validated parameters: `Y`, `Z`, `V`, `Tr`, `d`, and `r`.
#' @details
#' - The function validates and processes the input parameters:
#'   - \code{Y}, \code{Z}, \code{V}, and \code{Tr} default to \code{NA} if not provided.
#'   - \code{Z} and \code{V} must be matrices with the same number of columns (\code{2} for 2D or \code{3} for 3D).
#'   - \code{Tr} must have dimensions consistent with \code{Z} and \code{V}: 
#'     - For 2D: \code{Tr} must have 3 columns (triangles).
#'     - For 3D: \code{Tr} must have 4 columns (tetrahedra).
#'   - \code{d} and \code{r} are validated to ensure they meet the required conditions.
#'
#' @keywords internal
m <- function(..., Y = NULL, Z = NULL, V = NULL, Tr = NULL, d = NULL, r = NULL) {
  # Validate 'd': Must be a numeric value ≥ 1
  #if (missing(d) || is.null(d) || !is.numeric(d) || (d < 1)) {
  #  stop("Argument 'd' must be a numeric value greater than or equal to 1.")
  #}
  if (is.null(d)) {
    d.new <- NULL
  } else if (!is.null(d) && length(d) > 0) {
    if (!is.numeric(d) || (d < 1)) {
      stop("Argument 'd' must be a numeric value greater than or equal to 1.")
    } else {
      # Round 'd' and ensure it's an integer
      d.new <- round(d)
      if (!isTRUE(all.equal(d.new, d))) {
        stop("Argument 'd' must be an integer.")
      }
    }
  }
  
  # Round 'd' and ensure it's an integer
  # d.new <- round(d)
  # if (!isTRUE(all.equal(d.new, d))) {
  #   stop("Argument 'd' must be an integer.")
  # }
  
  # Validate 'r': Must be numeric
  if (missing(r) || is.null(r) || !is.numeric(r)) {
    stop("Argument 'r' must be a numeric value.")
  }
  
  # Round 'r' and ensure it's an integer
  r.new <- round(r)
  if (!isTRUE(all.equal(r.new, r))) {
    stop("Argument 'r' must be an integer.")
  }
  
  # Handle missing values for Z, V, Tr, and Y
  Z <- if (missing(Z) || is.null(Z)) NA else Z
  V <- if (missing(V) || is.null(V)) NA else V
  Tr <- if (missing(Tr) || is.null(Tr)) NA else Tr
  Y <- if (missing(Y) || is.null(Y)) NA else Y
  
  # Check if Z, V, and Tr are provided and validate their dimensions
  if (!anyNA(Z) && !anyNA(V) && !anyNA(Tr)) {
    # Ensure Z and V are matrices and have the same number of columns (2 or 3)
    if (!is.matrix(Z) || !is.matrix(V)) {
      stop("'Z' and 'V' must be matrices.")
    }
    
    if (ncol(Z) != ncol(V)) {
      stop("'Z' and 'V' must have the same number of columns (2 or 3).")
    }
    
    # Validate Tr's dimensions based on Z and V
    if (ncol(Z) == 2) {
      if (!is.matrix(Tr) || ncol(Tr) != 3) {
        stop("For 2D data: 'Z' and 'V' must each have 2 columns to describe points and vertices, and 'Tr' must have 3 columns to represent triangles.")
      }
    } else if (ncol(Z) == 3) {
      if (!is.matrix(Tr) || ncol(Tr) != 4) {
        stop("For 3D data: 'Z' and 'V' must each have 3 columns to describe points and vertices, and 'Tr' must have 4 columns to represent tetrahedra.")
      }
    } else {
      stop("'Z', 'V', and 'Tr' must have dimensions corresponding to either 2D or 3D space.")
    }
  }

  ret <- list(Y = Y, Z = Z, V = V, Tr = Tr, d = d.new, r = r.new)
  return(ret)
}


