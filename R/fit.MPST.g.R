#' Multivariate Penalized Spline over Triangulation Model Fitting using Global Method
#'
#' This function conducts the model fitting for multivariate penlized spline over triangulation using global method.
#'
#' @importFrom Matrix Matrix
#' 
#' @param Y The response variable observed over the domain.
#' \cr
#' @param Z The cooridinates of dimension \code{n} by two. Each row is the coordinates of a point.
#' \cr
#' @param V The \code{nV} by two matrix of vertices of a triangulation, where \code{nV} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 5, and usually \code{d} is greater than one. -1 represents piecewise constant.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param lambda The tuning parameter -- default is \eqn{10^(-6,-5.5,-5,\ldots,5,5.5,6)}.
#' \cr
#' @param Hmtx The indicator of whether the smoothness matrix \code{H} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param Kmtx The indicator of whether the energy matrix \code{K} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param QR The indicator of whether a QR decomposition need to be performed on the smoothness matrix -- default is \code{TRUE}.
#' \cr
#' @param TA The indicator of whether the area of the triangles need to be calculated -- default is \code{TRUE}.
#' \cr
#' @return A list of vectors and matrice, including:
#' \item{gamma.hat}{The estimated spline coefficients.}
#' \item{lamc}{The tuning parameter selected by Generalized Cross Validation (GCV).}
#' \item{B}{The spline basis function of dimension \code{n} by \code{nT}*\code{{(d+1)(d+2)/2}}, where \code{n} is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. The length of points means the length of ordering indices of observation points. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.}
#' \item{Ind.inside}{A vector contains the indexes of all the points which are inside the triangulation.}
#' \item{H}{The smoothness matrix.}
#' \item{Q2}{The Q2 matrix after QR decomposition of the smoothness matrix \code{H}.}
#' \item{K}{The thin-plate energy function.}
#' \item{tria.all}{The area of each triangle within the given triangulation.}
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
#' rmse = sqrt(mean((Y - mfit$Yhat)^2, na.rm = TRUE)); rmse
#' @export
#' 

fit.MPST.g <- function(Y, Z, V, Tr, d = 5, r = 1, lambda = 10^seq(-6, 6, by = 0.5)) {
  
  # 1. Preparation: basis generation, smoothness conditions, and penalty function
  n <- length(Y)
  
  B <- as.matrix(basis(V, Tr, d, r, Z)$B)
  
  if (d < 1 | r < 0 | nrow(Tr) <= 1) {
    warning("The degree of Bernstein polynomials d has be greater than zero, the smoothness parameter r has be nonnegative!")
    H <- NA; Q2 <- NA;
  } else {
    H <- as.matrix(smoothness(V, Tr, d, r))
    Q2 <- qrH(H)
  }
  
  if (d < 1) {
    warning("The degree of Bernstein polynomials d has be greater than zero.")
    K <- NA
  } else {
    K <- as.matrix(energy(V, Tr, d))
  }
  
  # tria.all <- TArea(V, Tr, Z)
  
  # 2. Model fitting
  nq <- ncol(Q2)
  
  W <- as.matrix(B %*% Q2)
  WW <- crossprod(W, W)
  rhs <- crossprod(W, Y)
  D <- crossprod(t(crossprod(Q2, as.matrix(K))), Q2)
  D <- as.matrix(D)
  
  flag <- (rankMatrix(WW) < nq)
  if(!flag){
    Ainv <- chol(WW, pivot = TRUE)
    A <- solve(t(Ainv))
    ADA <- A %*% D %*% t(A)
    eigs <- eigen(ADA)
    Cval <- eigs$values
  }
  
  nl <- length(lambda)
  theta.all <- c(); gamma.all <- c(); beta.all <- c()
  res.all <- c(); sse.all <- c(); df.all <- c(); gcv.all <- c(); bic.all <- c()
  for(il in 1:nl){
    Lam <- lambda[il]
    Dlam <- Lam * D
    lhs <- WW + Dlam
    # chol2inv(chol()) < solve() < qr.solve(); 
    lhs.inv <- chol2inv(chol(lhs));
    # crossprod() < %*%;
    theta <- crossprod(t(lhs.inv), rhs)
    theta.all <- cbind(theta.all, theta)
    gamma <- crossprod(t(Q2), theta)
    gamma.all <- cbind(gamma.all, gamma)
    beta <- crossprod(t(W), theta)
    beta.all <- cbind(beta.all, beta)
    res <- Y - beta
    res.all <- cbind(res.all, res)
    sse <- sum(res^2)
    sse.all <- c(sse.all, sse)
    if (!flag) {
      df <- sum(1 / (1 + Cval * Lam))
    } else {
      Hmtx <- crossprod(t(crossprod(t(W), lhs.inv)), t(W))
      df <- sum(diag(Hmtx))
    }
    df.all <- c(df.all, df)
    gcv <- n * sse / (n - df)^2
    gcv.all <- c(gcv.all, gcv)
    bic <- log(sse / n) + df * log(n) / n
    bic.all <- c(bic.all, bic)
  }
  j <- which.min(gcv.all)
  
  lambdac <- lambda[j]
  theta <- theta.all[, j]
  gamma <- gamma.all[, j]
  beta <- beta.all[, j]
  df <- df.all[j]
  sse <- sse.all[j]
  gcv <- gcv.all[j]
  bic <- bic.all[j]
  Yhat = beta
  res = Y - beta
  mse = mean((Y - Yhat)^2, na.rm = TRUE)
  
  mfit <- list(gamma.hat = gamma, 
               gamma.star = NULL,
               B = B, 
               B.star = NULL,
               lamc = lambdac, 
               sse = sse,
               mse = mse,
               gcv = gcv,
               bic = bic, 
               edf = df, 
               V = V,
               Tr = Tr)
  
  return(mfit)
}
