#' Make predictions from a fitted BPST object.
#'
#' This function is used to make predictions of a fitted BPST object.
#'
#' @importFrom Matrix Matrix
#' @importFrom pracma isempty
#' @importFrom Rcpp evalCpp
#'
#' @param mfit Fitted ``MPST" object.
#' \cr
#' @param Znew The cooridinates of new locationsfor prediction -- default is the observed coordinates, \code{Z}.
#' \cr
#' @return A vector of predicted values is returned.
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia and Li Wang from the Iowa State University.
#'
#' @examples
#' data(VT.square) 
#' d = 5; r = 1; 
#' func = 1; sigma = 1;
#' n = 2000;
#' Z = matrix(runif(2*n, 0, 1), nrow = n, ncol = 2)
#' sam = dataGenerator2D(Z, V, Tr, func, sigma)
#' Y = as.vector(sam$Y); Z = as.matrix(sam$Z);
#' mfit = fit.MPST(Y, Z, V, Tr, d, r)
#' # Grid points for evaluation
#' n1.grid = 51; n2.grid = 51; n.grid = n1.grid * n2.grid;
#' u.grid = seq(0, 1, length.out = n1.grid)
#' v.grid = seq(0, 1, length.out = n2.grid)
#' uu.grid = rep(u.grid, each = n2.grid)
#' vv.grid = rep(v.grid, times = n1.grid)
#' Z.grid = as.matrix(cbind(uu.grid, vv.grid))
#' pop = dataGenerator2D(Z.grid, V, Tr, func, sigma)
#' Y.grid = pop$Y; mu.grid = pop$mu;
#' mpred = pred.MPST(mfit, Z.grid)
#' rmspe = sqrt(mean((Y.grid - mpred$Ypred)^2, na.rm = TRUE)); rmspe
#'
#' @export
#'
pred.MPST <- function(mfit, Znew = NULL){
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
