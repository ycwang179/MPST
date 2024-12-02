#' Predict for MPST Models
#'
#' @description Predicts from a fitted MPST model using global ("G") or distributed ("D") methods.
#' @param formula A formula specifying the model to predict, e.g., `y ~ m(Z, V, Tr, d = 2, r = 1)`.
#' @param lambda A numeric vector specifying the smoothing parameters.
#' @param method A character string specifying the method: "G" (Global) or "D" (Distributed).
#' @param P.func An integer indicating the parallelization function: 1 (`mclapply`) or 2 (`parLapply`).
#' @param data A list containing required data for the model fitting.
#' @param data.pred A list containing prediction grid data, including `Z.grid`, `Y.grid`, and `mu.grid`.
#' @param debug Logical, if TRUE, prints additional debugging information.
#' @return An object of class "MPST" containing prediction results.
#' @export

predict.MPST <- function(formula, lambda, method, P.func = NULL, data = list(), data.pred = list()) {
  # Validate 'method'
  if (missing(method) || !(method %in% c("G", "D"))) {
    stop("Invalid 'method'. Please specify 'G' for Global or 'D' for Distributed learning.")
  }
  
  # Validate 'lambda'
  if (missing(lambda) || !is.numeric(lambda)) {
    stop("Invalid 'lambda'. Please provide a numeric vector of smoothing parameters.")
  }
  
  # Set default for 'P.func'
  if (is.null(P.func)) {
    P.func <- 2 # Default to parLapply
  } else if (!is.numeric(P.func) || !(P.func %in% c(1, 2))) {
    stop("Invalid 'P.func'. Use 1 for 'mclapply' or 2 for 'parLapply'.")
  }
  
  # Check required data components
  required_components <- c("Y", "Z", "V", "Tr")
  if (!all(required_components %in% names(data))) {
    stop("'data' must contain the following components: 'Y', 'Z', 'V', and 'Tr'.")
  }
  
  # Extract prediction data
  Z.grid <- data.pred$Z.grid
  mu.grid <- data.pred$mu.grid
  
  if (is.null(Z.grid)) stop("Missing required prediction grid: 'Z.grid'.")
  
  # Extract model data
  mpst.p <- data
  mpst.p$lambda <- lambda
  mpst.p$P.func <- P.func
  mpst.p$formula <- formula
  mpst.p$method <- method
  
  # Parse formula
  interp <- interpret.mpst(formula)
  mpst.p$d <- interp$d
  mpst.p$r <- interp$r

  # Fit and predict
  mfit <- fit.mpst(mpst.p$Y, mpst.p$Z, mpst.p$V, mpst.p$Tr, mpst.p$d, mpst.p$r, mpst.p$lambda, nl = 1, method = mpst.p$method, P.func = mpst.p$P.func)
  mpred <- pred.mpst(mfit, Z.grid)
  
  # Compute MISE
  mpred$mise <- mean((mu.grid - mpred$Ypred)^2, na.rm = TRUE)
  
  # Assign additional attributes
  mpred$P.func <- P.func
  mpred$method <- method
  mpred$formula <- formula
  mpred$func <- "predict"
  
  class(mpred) <- "MPST"
  
  return(mpred)
}
