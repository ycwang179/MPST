#' Model Fitting using Multiivariate Penalized Spline over Triangulation
#'
#' This function conducts the model fitting via multivariate penlized spline over triangulation.
#'
#' @importFrom Matrix Matrix
#' @importFrom MatrixExtra t_shallow t_deep
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

fit.MPST <- function(Y, Z, V, Tr, d = 5, r = 1, lambda = 10^seq(-6, 6, by = 0.5), nl = 1, method = "G") {
  
  this.call <- match.call()
  
  n <- length(Y)
  
  inVT.list = inVT(V, Tr, Z)
  ind.inside = which(inVT.list$ind.inside == 1)
  ind.nna <- (1:n)[!is.na(Y)]
  ind = intersect(ind.inside, ind.nna)
  Zi <- Z[ind, ]
  Yi <- Y[ind]
  ni = length(Yi)
  
  nd = ncol(Tr)
  if (nd == 3) {
    #nq = choose(d + 2, 2)
    i.max <- 4
  } else if (nd == 4) {
    #nq = choose(d + 3, 3)
    i.max <- 9
  }
  
  if (method == "G") {
    N.cores <- 1
    all.info <- list()
    best.list <- list()
    
    if ((!hasArg(d)) || is.null(d) || (d < 1)) {
      
      for (i in (1 : i.max)) {
        d <- i + 1
        mfit.temp <- fit.MPST.g(Yi, Zi, V, Tr, d, r, lambda)
        all.info[[length(all.info) + 1]] <- list(d = d, mfit.temp)
      }
      
      gcv.d.all <- sapply(all.info, function(x) x[[2]]$gcv)
      index.d.gcv <- which.min(gcv.d.all)
      best.list <- all.info[[index.d.gcv]]
      
    } else {
      
      d <- d
      mfit.temp <- fit.MPST.g(Yi, Zi, V, Tr, d, r, lambda)
      best.list <- list(d = d, mfit.temp)
      
    }
    
    mfit <- unlist(best.list, recursive = FALSE)
    
    all.info <- list()
    best.list <- list()
    
    d = mfit$d
    gamma.hat = mfit$gamma.hat
    gamma.star = NULL
    B = mfit$B
    B.star = mfit$B.star
    lambdac = mfit$lamc
    
  } else if (method == "D") {
    ns <- parallel::detectCores()
    N.cores <- ns
    #ns = 16
    
    if ((!hasArg(d)) || is.null(d) || (d < 1)) {
      d <- 5
    }
    
    if (nd == 3) {
      nq = choose(d + 2, 2)
      TV = as.matrix(tdata(V, Tr)$TV)
    } else if (nd == 4) {
      nq = choose(d + 3, 3)
      TV = as.matrix(thdata(V, Tr)$TV)
    }
    
    #load.all = worker.load(V = V, Tr = Tr, TV = TV, inVT.list = inVT.list, 
    #                       Y = Yi, Z = Zi, d = d, nl = nl, ns = ns)
    
    load.all = worker.load(V = V, Tr = Tr, TV = TV, inVT.list = inVT.list, 
                           Y = Yi, Z = Zi, d = d, nl = nl, ns = ns, P.func = P.func)
    
    if (P.func == 1) {
      mfit.all <- parallel::mclapply(1:nrow(Tr), FUN = fit.MPST.d, mc.cores = ns,
                                     Y = Yi, Z = Zi, V = V, Tr = Tr, d = d, r = r, 
                                     lambda = lambda, nl = nl, load.all = load.all)
    } else if (P.func == 2) {
      cl <- parallel::makeCluster(ns)
      
      # Load the necessary packages on each worker node
      parallel::clusterEvalQ(cl, {
        library(pracma)
        library(Matrix)
      })
      
      # Export custom functions and variables to the cluster
      parallel::clusterExport(cl, varlist = c("fit.MPST.d", "n","Yi", "Zi", "V", "Tr", "d", "r", "lambda", "nl", "load.all", "mtxcbind"), envir = environment())
      
      # Replace mclapply with parLapply
      mfit.all <- parallel::parLapply(cl, 1:nrow(Tr), function(iT) {
        fit.MPST.d(iT, Y = Yi, Z = Zi, V = V, Tr = Tr, d = d, r = r, 
                   lambda = lambda, nl = nl, load.all = load.all)
      })
      
      # Stop the parallel cluster
      parallel::stopCluster(cl)
      
    }
    
    # mfit.all = vector(mode = "list", length = nrow(Tr))
    # for (iT in 1:nrow(Tr)) {
    #   mfit.all[[iT]] = fit.MPST.d(iT, Y, Z, V, Tr, d, r, lambda, nl, load.all)
    # }
    
    # combine results from workers
    for (argument in names(mfit.all[[1]])) {
      assign(argument, lapply(mfit.all, "[[", argument))
    }
    
    # re-order gamma
    gamma.star <- unlist(gamma.star)
    
    # linear coefficients gamma for bivariate function
    if (nd == 3) {
      H = as.matrix(smoothness(V, Tr, d, r))
    } else if (nd == 4) {
      H = as.matrix(smoothness3D(V, Tr, d, r))
    }
    
    a = H %*% gamma.star
    HH = crossprod(t(H)); nH = nrow(HH)
    b = chol2inv(chol(HH + 1e-12 * diag(nH))) %*% a
    gamma.hat = gamma.star - t(H) %*% b
    
    B = matrix(0, nrow = ni, ncol = nq * nrow(Tr))
    for (iT in 1:nrow(Tr)) {
      B[inVT.list$ind.T == iT, (iT - 1) * nq + 1:nq] = B.star[[iT]]
    }
  }
  
  Y.hat = rep(NA, n)
  Y.hat[ind] = B %*% gamma.hat
  res = Y - Y.hat
  sse = sum(res^2, na.rm = TRUE)
  mse = mean(res^2, na.rm = TRUE)
  
  mfit <- list(gamma.hat = gamma.hat, 
               gamma.star = gamma.star, 
               B = B,
               B.star = B.star,
               Y.hat = Y.hat, 
               lamc = lambdac, 
               res = res,
               sse = sse,
               mse = mse,
               V = V,
               Tr = Tr,
               d = d,
               r = r,
               Y = Y,
               Z = Z,
               N.cores = N.cores)
  
  mfit$call <- this.call;
  class(mfit) <- "MPST"
  
  return(mfit)
}
