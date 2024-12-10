#' Predict for MPST Models
#'
#' @description This method generates predictions from a fitted MPST model using global (`"G"`) 
#' or distributed (`"D"`) learning methods. It allows for flexible prediction grids and computes 
#' the mean integrated squared error (MISE) if a reference function (`mu.grid`) is provided.
#'
#' @rdname predict
#' @method predict MPST
#' @param formula A formula specifying the model, e.g., `y ~ m(Z, V, Tr, d, r)`. 
#' - `Y`: The response variable observed over the domain.
#' - `Z`: Matrix of observation coordinates (\code{n} by \code{k}). Rows represent points in 
#'   2D or 3D space (\code{k = 2} or \code{k = 3}). \( k \) is the dimension of the observed 
#'   points, where \( k = 2 \) for 2D and \( k = 3 \) for 3D.
#' - `V`: Matrix of vertices (\code{nV} by \code{k}). Rows represent coordinates of vertices 
#'   in the triangulation.
#' - `Tr`: Triangulation matrix (\code{nT} by \code{k+1}). Rows represent vertex indices:
#'   - For 2D: Rows have three indices for triangles.
#'   - For 3D: Rows have four indices for tetrahedra.
#' - `d`: Degree of piecewise polynomials (default: \code{5}). \code{-1} represents piecewise constants.
#' - `r`: Smoothness parameter (default: \code{1}, where \code{0 <= r < d}).
#'
#' @param lambda The tuning parameter. If not specified, defaults to \eqn{10^(-6,-5.5,-5,\ldots,5,5.5,6)}.
#' @param method A character string specifying the learning method. If not specified, defaults to `"G"` (Global learning).
#' - `"G"`: Global learning.
#' - `"D"`: Distributed learning.
#' @param P.func An integer specifying the parallelization method for distributed learning. Defaults to \code{2}:
#' - `1`: Use `mclapply`.
#' - `2`: Use `parLapply`.
#' @param data (Optional) A list containing the following components:
#' - `Y`: The response variable observed over the domain.
#' - `Z`: Matrix of observation coordinates.
#' - `V`: Matrix of triangulation vertices.
#' - `Tr`: Triangulation matrix.
#' 
#' @param data.pred A list containing prediction-related data:
#' - `Z.grid`: The prediction grid coordinates (required).
#' - `mu.grid`: (Optional) True mean values for computing mean integrated squared error (MISE).
#'
#' @return An object of class `"MPST"` containing the following components:
#' - `Ypred`: Predicted values on the specified grid.
#' - `mise`: Mean integrated squared error (computed if `mu.grid` is provided).
#' - `method`: Learning method used for prediction.
#' - `formula`: The formula used for fitting and prediction.
#' 
#' @examples
#' \dontrun{
#' # Example using a fitted model and prediction grid
#' data_list <- list(Y = y, Z = Z, V = V, Tr = Tr)
#' data_pred <- list(Z.grid = Z_new, mu.grid = mu_new)
#' predictions <- predict.MPST(y ~ m(Z, V, Tr, d = 2, r = 1), data = data_list, data.pred = data_pred)
#' 
#' # View predictions
#' print(predictions$Ypred)
#' # If MISE was computed
#' print(predictions$mise)
#' }
#' @export
predict.MPST <- function(formula, lambda = NULL, method = NULL, P.func = NULL, data = list(), data.pred = list()) {

  `%||%` <- function(a, b) if (!is.null(a)) a else b
    
  # Check if the parameter formula is provided
  if (missing(formula)) {
    stop("'formula' is required. Please specify a formula (e.g., y ~ m(Z, V, Tr, d, r)).")
  }

  # Set default parameters and check validity
  method <- method %||% "G"  # Default to Global Learning
  if (!(method %in% c("G", "D"))) {
    stop("Invalid 'method'. Use 'G' for Global or 'D' for Distributed learning.")
  }
  
  lambda <- lambda %||% 10^seq(-6, 6, by = 0.5) # Default range for lambda
  if (!is.numeric(lambda)) {
    stop("Invalid 'lambda'. Please provide a numeric vector of smoothing parameters.")
  }
  
  P.func <- P.func %||% 2 # Default to parLapply
  if (!is.numeric(P.func) || !(P.func %in% c(1, 2))) {
    stop("Invalid 'P.func'. Use 1 for 'mclapply' or 2 for 'parLapply'.")
  }

  # Check if the data contains the required components
  required_components <- c("Y", "Z", "V", "Tr")
  missing_components <- setdiff(required_components, names(data))
  if (length(missing_components) > 0) {
    stop(paste0("'data' must contain the following components: ", 
                paste(required_components, collapse = ", "), 
                ". Missing components: ", 
                paste(missing_components, collapse = ", ")))
  }
  
  # Extract prediction data
  if (!("Z.grid" %in% names(data.pred))) {
    stop("Missing required prediction grid: 'Z.grid'.")
  }
  Z.grid <- data.pred$Z.grid
  mu.grid <- data.pred$mu.grid
  
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
  mfit <- fit.mpst.internal(mpst.p$Y, mpst.p$Z, mpst.p$V, mpst.p$Tr, mpst.p$d, mpst.p$r, mpst.p$lambda, nl = 1, method = mpst.p$method, P.func = mpst.p$P.func)
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
