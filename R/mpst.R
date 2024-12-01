#' Print Method for MPST Objects
#'
#' @description Provides a custom print method for objects of class "MPST".
#' @param x An object of class "MPST".
#' @param digits Number of significant digits to display (default: 4).
#' @param max Maximum number of elements to print for long vectors (default: 10).
#' @return Invisibly returns the input object `x`.
#' @examples
#' \dontrun{
#' print.MPST(mpst_object)
#' }
#' @export
print.MPST <- function(x, digits = 4, max = 10) {
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
    cat("\nMSE:", round(x$mse, digits), "\n")  
    cat("\ngamma_hat:\n")
    print(head(x$gamma.hat, max))
    cat("\ny_hat:\n")
    print(head(x$Y.hat, max))
    
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
    cat("\nMISE:", round(x$mise, digits), "\n")  
    cat("\ny_predict:\n")
    print(head(x$Ypred, max))
    
  } else if (x$func == "basis") {
    # handle basis Situations
    if (x$method == "G") {
      cat("\nMethod: Global\n")
      cat("\nbasis:\n")
      print(head(x$B, max))  
    } else if (x$method == "D") {
      cat("\nMethod: Distributed Learning\n")
      for (item in 1:length(x$basis.all)) {
        cat("\nbasis for triangle", item, ":\n")
        print(head(x$basis.all[[item]]$B.star, max))  
      }
    }
  }
  
  invisible(x)
}

#' Interpret an MPST Formula
#'
#' @description Parses and interprets an MPST formula for use in model fitting.
#' @param mpstf An MPST formula object.
#' @param extra.special Additional special terms to recognize (default: NULL).
#' @return A list containing parsed terms and the response variable.
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

#' Generate MPST Parameters
#'
#' @description Internal function to create a list of parameters for MPST based on input arguments.
#' This function is not intended for direct use by users.
#' @param ... Additional arguments (not used).
#' @param V Matrix of vertices (optional).
#' @param Tr Matrix of triangulations (optional).
#' @param d Degree of the spline (must be an integer ≥ 1).
#' @param r Smoothness parameter (must be an integer).
#' @return A list containing parameters `d` and `r`.
#' @keywords internal
m <- function(..., V = NULL, Tr = NULL, d, r) {
  if (missing(d) || is.null(d) || d < 1) {
    d.new <- NULL
  } else {
    d.new <- round(d)
    if (!isTRUE(all.equal(d.new, d))) {
      stop("Argument 'd' should be an integer.")
    }
  }
  
  r.new <- round(r)
  if (!isTRUE(all.equal(r.new, r))) {
    stop("Argument 'r' should be an integer.")
  }
  
  ret <- list(d = d.new, r = r.new)
  return(ret)
}


