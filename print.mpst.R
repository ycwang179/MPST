#' Print Method for MPST Class
#'
#' This function provides a customized print method for objects of class \code{MPST}. 
#' Depending on the object type (\code{fit}, \code{predict}, or \code{basis}), 
#' it prints relevant information, such as formulas, methods, fitted values, predictions, and basis functions.
#'
#' @param x An object of class \code{MPST}.
#' @param digits Integer. The number of decimal places to round numeric values. Default is 4.
#' @param max Integer. The maximum number of elements to display in printed vectors or matrices. Default is 10.
#' 
#' @return The input object \code{x}, invisibly.
#' @examples
#' # Assume `obj` is an object of class MPST
#' # print.MPST(obj, digits = 3, max = 5)
#' 
#' # Use the S3 generic `print`:
#' # print(obj)
#' 
#' @export
#' @method print MPST
print.MPST <- function(x, digits = 4, max = 10) {
  if (!inherits(x, "MPST")) {
    stop("Object must be of class 'MPST'")
  }
  
  cat("\nFormula: ")
  print(x$formula)
  
  if (x$func == "fit") {
    # Handle fit Situations
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
    # Handle predict Situations
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
    # Handle basis Situations
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
