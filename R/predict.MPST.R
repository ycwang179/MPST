#' Predict for MPST Models
#'
#' @description Predicts from a fitted MPST model using global ("G") or distributed ("D") methods.
#'
#' @param formula A formula specifying the model, e.g., `y ~ m(Z, V, Tr, d, r)`. 
#' - `Y`: The response variable observed over the domain.
#' - `Z`: Matrix of observation coordinates (n by k). Rows represent points in 
#'   2D or 3D space (k = 2 or k = 3).
#' - `V`: Matrix of vertices (nV by k). Rows represent coordinates of vertices 
#'   in the triangulation.
#' - `Tr`: Triangulation matrix (nT by k+1). Rows represent vertex indices:
#'   - For 2D: Rows have three indices for triangles.
#'   - For 3D: Rows have four indices for tetrahedra.
#' - `d`: Degree of piecewise polynomials.
#' - `r`: Smoothness parameter.
#'
#' @param lambda The tuning parameter -- default is \eqn{10^(-6,-5.5,-5,\ldots,5,5.5,6)}.
#' @param method A character string specifying the learning method:
#' - `"G"`: Global learning.
#' - `"D"`: Distributed learning.
#' @param P.func An integer specifying the parallelization method for distributed learning (default: 2):
#' - `1`: Use `mclapply`.
#' - `2`: Use `parLapply`.
#' @param data A list containing the following components:
#' - `Y`: The response variable observed over the domain.
#' - `Z`: Matrix of observation coordinates.
#' - `V`: Matrix of triangulation vertices.
#' - `Tr`: Triangulation matrix.
#' @param data.pred A list containing prediction grid data:
#' - `Z.grid`: The prediction grid coordinates.
#' - `mu.grid`: (Optional) True mean values for computing mean integrated squared error (MISE).
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
  
  # Compute MISE if 'mu.grid' is provided
  if (!is.null(mu.grid)) {
    mpred$mise <- mean((mu.grid - mpred$Ypred)^2, na.rm = TRUE)
  }
  
  # Assign additional attributes
  mpred$P.func <- P.func
  mpred$method <- method
  mpred$formula <- formula
  mpred$func <- "predict"
  
  class(mpred) <- "MPST"
  
  return(mpred)
}

#' Predict Values for New Data
#'
#' @description Internal function for making predictions using a fitted MPST model.
#'
#' @importFrom Matrix Matrix
#' @importFrom pracma isempty
#' @importFrom Rcpp evalCpp
#'
#' @param mfit A fitted MPST model object.
#' @param Znew A matrix of new coordinates for prediction.
#' @return A list containing:
#' - `Ypred`: Predicted values.
#' - `ind.inside`: Indices of the points inside the domain.
#' - `d`: Degree of the spline.
#' - `N.cores`: Number of cores used for computation.
#' @keywords internal
pred.mpst <- function(mfit, Znew = NULL){
  if(identical(Znew, mfit$Z) | isempty(Znew)){
    Ypred <- mfit$beta.hat
    ind.inside <- mfit$ind.inside
  }else{
    nd = ncol(mfit$Tr)
    if (nd == 3) {
      B.all <- basis2D.d(mfit$V, mfit$Tr, mfit$d, mfit$r, Znew)
      Bnew = B.all$B
      ind.inside <- B.all$ind.inside
    } else if (nd == 4) {
      B.all <- basis3D.d(mfit$V, mfit$Tr, mfit$d, mfit$r, Znew)
      Bnew = B.all$B
      ind.inside <- B.all$ind.inside
    }
    Ypred <- rep(NA, nrow(Znew))
    Ypred[ind.inside] <- Bnew %*% mfit$gamma.hat
  }
  mpred = list(Ypred = Ypred, 
               ind.inside = ind.inside,
               d = mfit$d,
               N.cores = mfit$N.cores)
  return(mpred)
}
