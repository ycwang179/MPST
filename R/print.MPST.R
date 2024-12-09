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
