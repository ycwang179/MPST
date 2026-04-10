#' Predict for MPST Models
#'
#' @description `predict.MPST()` generates predictions for a Multivariate Penalized Spline
#' over Triangulation (MPST) model using global (`"G"`) or distributed (`"D"`) learning.
#' It allows prediction on user-supplied locations and can compute the mean integrated squared
#' error (MISE) when a reference mean function is provided.
#'
#' @rdname predict
#' @method predict MPST
#' @param formula A formula specifying the model, e.g., `Y ~ m(Z, V, Tr, d, r)`.
#' - `Y`: The response variable observed over the domain.
#' - `Z`: Matrix of observation coordinates (\code{n} by \code{k}). Rows represent points in
#'   2D or 3D space (\code{k = 2} or \code{k = 3}).
#' - `V`: Matrix of vertices (\code{nV} by \code{k}). Rows represent coordinates of vertices
#'   in the triangulation.
#' - `Tr`: Triangulation matrix (\code{nT} by \code{k+1}). Rows represent vertex indices:
#'   - For 2D: Rows have three indices for triangles.
#'   - For 3D: Rows have four indices for tetrahedra.
#' - `d`: Degree of piecewise polynomials. If `d = NULL`, the degree is selected automatically
#'   according to the learning method. Under global learning (`method = "G"`), GCV is used to
#'   select the degree from \code{2:5} for 2D triangulations and from \code{2:9} for 3D
#'   triangulations. Under distributed learning (`method = "D"`), the default degree is
#'   \code{5}. `-1` represents piecewise constants.
#' - `r`: Smoothness parameter (default: \code{1}, where \code{0 <= r < d}).
#'
#' @param lambda A numeric vector of tuning parameters for regularization. Defaults to
#' \eqn{10^(-6,-5.5,-5,\ldots,5,5.5,6)}.
#' @param method A character string specifying the learning method. If not specified,
#' defaults to `"G"` (global learning).
#' - `"G"`: Global learning.
#' - `"D"`: Distributed learning.
#' @param P.func An integer specifying the parallelization method for distributed learning.
#' Defaults to \code{2}:
#' - `1`: Use `mclapply`.
#' - `2`: Use `parLapply`.
#' @param data (Optional) A list containing the following components:
#' - `Y`: The response variable observed over the domain.
#' - `Z`: Matrix of observation coordinates.
#' - `V`: Matrix of triangulation vertices.
#' - `Tr`: Triangulation matrix.
#' - `d`: (Optional) Degree of piecewise polynomials. If not supplied in either `formula`
#'   or `data`, or if `d = NULL`, the function uses method-specific default behavior:
#'   under global learning, GCV selects `d` from \code{2:5} in 2D and \code{2:9} in 3D;
#'   under distributed learning, `d = 5` is used.
#' - `r`: Smoothness parameter.
#'
#' @param data.pred A list containing prediction-related data:
#' - `Z.grid`: The prediction coordinates (required).
#' - `mu.grid`: (Optional) True mean values used to compute mean integrated squared error (MISE).
#'
#' @return An object of class `"MPST"` containing:
#' - `Ypred`: Predicted values at the specified grid locations.
#' - `ind.inside`: Indices of prediction points lying inside the domain.
#' - `mise`: Mean integrated squared error, if `mu.grid` is provided.
#' - `method`: The learning method used for prediction.
#' - `formula`: The formula used for fitting and prediction.
#'
#' @details
#' This function first fits an MPST model using the supplied `formula`, `data`, and tuning
#' parameters, and then generates predictions on `data.pred$Z.grid`.
#'
#' The degree parameter `d` is optional. When `d = NULL`, the function uses the same automatic
#' degree selection rule as `fit.MPST()`: under global learning, GCV selects the degree from
#' \code{2:5} in 2D and \code{2:9} in 3D; under distributed learning, the degree defaults to
#' \code{5}.
#'
#' If a required component other than `d` is missing in both `formula` and `data`, the function
#' raises an error. The prediction grid `Z.grid` must be provided through `data.pred`.
#'
#' @export
predict.MPST <- function(formula, lambda = NULL, method = NULL, P.func = NULL, data = list(), data.pred = list()) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b
    
  if (missing(formula)) {
    stop("'formula' is required. Please specify a formula (e.g., y ~ m(Z, V, Tr, d, r)).")
  }

  method <- method %||% "G"
  if (!(method %in% c("G", "D"))) {
    stop("Invalid 'method'. Use 'G' for Global or 'D' for Distributed learning.")
  }
  
  lambda <- lambda %||% 10^seq(-6, 6, by = 0.5)
  if (!is.numeric(lambda)) {
    stop("Invalid 'lambda'. Please provide a numeric vector of smoothing parameters.")
  }
  
  P.func <- P.func %||% 2
  if (!is.numeric(P.func) || !(P.func %in% c(1, 2))) {
    stop("Invalid 'P.func'. Use 1 for 'mclapply' or 2 for 'parLapply'.")
  }
  
  interp <- interpret.mpst(formula)
  
  Y <- if (!is.null(interp$Y) && !all(is.na(interp$Y))) interp$Y else data$Y
  Z <- if (!is.null(interp$Z) && !all(is.na(interp$Z))) interp$Z else data$Z
  V <- if (!is.null(interp$V) && !all(is.na(interp$V))) interp$V else data$V
  Tr <- if (!is.null(interp$Tr) && !all(is.na(interp$Tr))) interp$Tr else data$Tr
  d <- if (!is.null(interp$d) && !all(is.na(interp$d))) interp$d else data$d
  r <- if (!is.null(interp$r) && !all(is.na(interp$r))) interp$r else data$r

  if (is.null(Y) || is.null(Z) || is.null(V) || is.null(Tr) || is.null(r)) {
    stop("Both 'formula' and 'data' must provide the components: 'Y', 'Z', 'V', 'Tr', and 'r'.")
  }
  
  mpst.p <- list(
    Y = Y,
    Z = Z,
    V = V,
    Tr = Tr,
    d = d,
    r = r,
    lambda = lambda,
    P.func = P.func,
    method = method,
    formula = formula
  )
 
  if (!("Z.grid" %in% names(data.pred))) {
    stop("Missing required prediction grid: 'Z.grid'.")
  }
  Z.grid <- data.pred$Z.grid
  mu.grid <- data.pred$mu.grid
  
  mfit <- fit.mpst.internal(
    mpst.p$Y, mpst.p$Z, mpst.p$V, mpst.p$Tr,
    mpst.p$d, mpst.p$r, mpst.p$lambda,
    nl = 1, method = mpst.p$method, P.func = mpst.p$P.func
  )
  mpred <- pred.mpst(mfit, Z.grid)
  
  if (!is.null(mu.grid)) {
    mpred$mise <- mean((mu.grid - mpred$Ypred)^2, na.rm = TRUE)
  }
  
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
