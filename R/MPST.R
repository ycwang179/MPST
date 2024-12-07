# MPST.R - Core Functions and Generic Definitions

#' @export
fit <- function(formula, lambda = NULL, method = NULL, P.func = NULL, data = list()) {
  UseMethod("fit")
}

#' @export
predict <- function(formula, lambda = NULL, method = NULL, P.func = NULL, data = list(), data.pred = list()) {
  UseMethod("predict")
}

# Define the print.MPST function
#' Print Method for MPST Models
#'
#' @importFrom utils head
#'
#' @description Prints the summary of an MPST model.
#' @param x An MPST object.
#' @return No return value; prints to the console.
#' @export
print.MPST <- function(x, ...) {
  if (!inherits(x, "MPST")) {
    stop("Object must be of class 'MPST'")
  }
  
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
    print(head(x$gamma.hat, 10))
    cat("\ny_hat:\n")
    print(head(x$Y.hat, 10))
    
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
    print(head(x$Ypred, 10))
    
  } else if (x$func == "basis") {
    # handle basis Situations
    if (x$method == "G") {
      cat("\nMethod: Global\n")
      cat("\nbasis:\n")
      print(head(x$B, 10))  
    } else if (x$method == "D") {
      cat("\nMethod: Distributed Learning\n")
      for (item in 1:length(x$basis.all)) {
        cat("\nbasis for triangle", item, ":\n")
        print(head(x$basis.all[[item]]$B.star, 10))  
      }
    }
  }
  
  invisible(x)
}

# Define the interpret.mpst function
#' Interpret an MPST Formula
#'
#' @description Parses and interprets an MPST formula for use in model fitting.
#' @param mpstf An MPST formula object.
#' @param extra.special Additional special terms to recognize (default: NULL).
#' @return A list containing:
#'   - `terms`: Parsed formula terms.
#'   - `response`: Response variable (if present).
#' @keywords internal
interpret.mpst <- function(mpstf, extra.special = NULL) {
  mpst.tf <- terms.formula(mpstf, specials = c("m", extra.special)) 
  mpst.terms <- attr(mpst.tf, "term.labels") 
  
  if (attr(mpst.tf, "response") > 0) {  
    mpst.response <- as.character(attr(mpst.tf, "variables")[2])
  } else { 
    mpst.response <- NULL
  }
  mpst.parameters <- eval(parse(text = paste(mpst.terms[1])))
}

# Define the m function
#' Generate MPST Parameters
#'
#' @description Internal function to create a list of parameters for MPST.
#' @param ... Additional arguments (not used).
#' @param V Matrix of vertices (optional).
#' @param Tr Matrix of triangulations (optional).
#' @param d Degree of the spline (must be an integer ≥ 1).
#' @param r Smoothness parameter (must be an integer).
#' @return A list containing parameters `d` and `r`.
#' @keywords internal
m <- function(..., V = NULL, Tr = NULL, d, r) {
  # Validate d
  if (missing(d) || is.null(d) || !is.numeric(d) || (d < 1)) {
    stop("Argument 'd' must be a numeric value greater than or equal to 1.")
  }
  
  # Round d and ensure it's an integer
  d.new <- round(d)
  if (!isTRUE(all.equal(d.new, d))) {
    stop("Argument 'd' must be an integer. It was rounded, but this is not allowed.")
  }
  
  # Validate r
  if (missing(r) || is.null(r) || !is.numeric(r)) {
    stop("Argument 'r' must be a numeric value.")
  }
  
  # Round r and ensure it's an integer
  r.new <- round(r)
  if (!isTRUE(all.equal(r.new, r))) {
    stop("Argument 'r' must be an integer. It was rounded, but this is not allowed.")
  }
  
  # Return the validated list
  ret <- list(d = d.new, r = r.new)
  return(ret)
}


