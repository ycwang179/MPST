#' Fit a Multivariate Penalized Spline Model
#'
#' @description `fit.MPST()` fits a Multivariate Penalized Spline over Triangulation (MPST)
#' model using global (`"G"`) or distributed (`"D"`) learning methods. The function accepts
#' both a formula and a list of data, along with optional parameters for customization.
#'
#' @rdname fit
#' @method fit MPST
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
#' @param lambda The tuning parameter. If not specified, defaults to \eqn{10^(-6,-5.5,-5,...,5,5.5,6)}.
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
#' @return An object of class `"MPST"` with the following components:
#' - `gamma.hat`: Estimated spline coefficients from the fitted model.
#' - `Y.hat`: Predicted values based on the fitted model.
#' - `mse`: Mean squared error of the model (computed during fitting).
#' - `mise`: Mean integrated squared error.
#' - `method`: The learning method used ("G" or "D").
#' - `formula`: The formula provided during fitting.
#'
#' @export
fit.MPST <- function(formula, lambda = NULL, method = NULL, P.func = NULL, data = list()) {
  
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
  
  # Extract Y, Z, V, Tr, d, and r from formula or data.
  interp <- interpret.mpst(formula)
  Y <- interp$Y %||% data$Y
  Z <- interp$Z %||% data$Z
  V <- interp$V %||% data$V
  Tr <- interp$Tr %||% data$Tr
  d <- interp$d %||% data$d
  r <- interp$r %||% data$r
  
  # Check for the presence of required parameters.
  if (is.null(Y) || is.null(Z) || is.null(V) || is.null(Tr) || is.null(d) || is.null(r)) {
    stop("Both 'formula' and 'data' must provide the components: 'Y' 'Z', 'V', 'Tr', 'd', and 'r'.")
  }
  
  # Construct the parameter list.
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
  # Fit the model
  mfit <- fit.mpst.internal(mpst.p$Y,mpst.p$Z, mpst.p$V, mpst.p$Tr, mpst.p$d, mpst.p$r, mpst.p$lambda, nl = 1, mpst.p$method, P.func = mpst.p$P.func)
  
  mfit$P.func <- mpst.p$P.func
  mfit$method <- mpst.p$method
  mfit$formula <- mpst.p$formula
  mfit$func <- "fit"
  class(mfit) <- "MPST"
  return(mfit)
}

#' @name fit.mpst.internal
#' @title Internal Function for MPST Model Fitting
#' @description Internal function to fit a multivariate penalized spline over triangulation 
#' (MPST) using global or distributed learning.
#'
#' @importFrom Matrix Matrix
#' @importFrom MatrixExtra t_shallow t_deep
#'
#' @param Y The response variable vector.
#' @param Z A matrix of observation coordinates.
#' @param V A matrix of vertices of the triangulation.
#' @param Tr A triangulation matrix.
#' @param d The degree of piecewise polynomials.
#' @param r The smoothness parameter.
#' @param lambda The tuning parameter.
#' @param nl Number of lambda values to consider.
#' @param method Learning method: `"G"` for global or `"D"` for distributed.
#' @param P.func Parallelization method for distributed learning.
#' @return A list containing model fit components.
#' @keywords internal
fit.mpst.internal <- function(Y, Z, V, Tr, d = NULL, r = 1, lambda, nl = 1, method, P.func) {
  
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
    i.max <- 8
  }
  
  if (method == "G") {
    N.cores <- 1
    all.info <- list()
    best.list <- list()
    
    if ((!hasArg(d)) || is.null(d) || (d < 1)) {
      
      for (i in (1 : i.max)) {
        d <- i + 1
        mfit.temp <- fit.mpst.g(Yi, Zi, V, Tr, d, r, lambda)
        all.info[[length(all.info) + 1]] <- list(d = d, mfit.temp)
      }
      
      gcv.d.all <- sapply(all.info, function(x) x[[2]]$gcv)
      index.d.gcv <- which.min(gcv.d.all)
      best.list <- all.info[[index.d.gcv]]
      
    } else {
      
      d <- d
      mfit.temp <- fit.mpst.g(Yi, Zi, V, Tr, d, r, lambda)
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
      mfit.all <- parallel::mclapply(1:nrow(Tr), FUN = fit.mpst.d, mc.cores = ns,
                                     Y = Yi, Z = Zi, V = V, Tr = Tr, d = d, r = r, 
                                     lambda = lambda, nl = nl, load.all = load.all)
    } else if (P.func == 2) {
      cl <- parallel::makeCluster(ns)
      
      # Load the necessary packages on each worker node
      parallel::clusterEvalQ(cl, {
        library(pracma)
        library(Matrix)
        library(MatrixExtra) 
      })
      
      # Export custom functions and variables to the cluster
      parallel::clusterExport(cl, varlist = c("fit.mpst.d", "n","Yi", "Zi", "V", "Tr", "d", "r", "lambda", "nl", "load.all", "mtxcbind"), envir = environment())
      
      mfit.all <- parallel::parLapply(cl, 1:nrow(Tr), function(iT) {
        fit.mpst.d(iT, Y = Yi, Z = Zi, V = V, Tr = Tr, d = d, r = r, 
                   lambda = lambda, nl = nl, load.all = load.all)
      })
      
      # Stop the parallel cluster
      parallel::stopCluster(cl)
      
    }
    
    # mfit.all = vector(mode = "list", length = nrow(Tr))
    # for (iT in 1:nrow(Tr)) {
    #   mfit.all[[iT]] = fit.mpst.d(iT, Y, Z, V, Tr, d, r, lambda, nl, load.all)
    # }
    
    # combine results from workers
    for (argument in names(mfit.all[[1]])) {
      assign(argument, lapply(mfit.all, "[[", argument))
    }
    
    # re-order gamma
    gamma.star <- unlist(gamma.star)
    
    # linear coefficients gamma for bivariate function
    if (nd == 3) {
      H = as.matrix(smoothness2D(V, Tr, d, r))
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

#' Fit an MPST Model Using Global Learning
#'
#' @description Internal function to fit a multivariate penalized spline over 
#' triangulation (MPST) using global learning.
#'
#' @importFrom Matrix Matrix rankMatrix
#'
#' @param Y The response variable vector.
#' @param Z A matrix of observation coordinates.
#' @param V A matrix of vertices of the triangulation.
#' @param Tr A triangulation matrix.
#' @param d The degree of piecewise polynomials.
#' @param r The smoothness parameter.
#' @param lambda The tuning parameter.
#' @return A list containing global model fit components.
#' @keywords internal
fit.mpst.g <- function(Y, Z, V, Tr, d = NULL, r = 1, lambda) {
  
  # 1. Preparation: basis generation, smoothness conditions, and penalty function
  n <- length(Y)
  nd = ncol(Tr)
  if (d < 1 | r < 0 | nrow(Tr) <= 1) {
    warning("The degree of Bernstein polynomials d has be greater than zero, the smoothness parameter r has be nonnegative!")
    B <- NA
    H <- NA; Q2 <- NA;
    K <- NA
  } else {
    if (nd == 3) {
      B <- basis2D.d(V, Tr, d, r, Z)$B
      H <- as.matrix(smoothness2D(V, Tr, d, r))
      K <- energy(V, Tr, d)
    } else if (nd == 4) {
      B <- basis3D.d(V, Tr, d, r, Z)$B
      H <- as.matrix(smoothness3D(V, Tr, d, r))
      K <- energy3D(V, Tr, d)
    }
    Q2 <- qrH(H)
  }
  # tria.all <- TArea(V, Tr, Z)
  
  # 2. Model fitting
  nq <- ncol(Q2)
  
  W <- as.matrix(B %*% Q2)
  WW <- crossprod(W, W)
  rhs <- crossprod(W, Y)
  D <- crossprod(t(crossprod(Q2, as.matrix(K))), Q2)
  D <- as.matrix(D)
  
  flag <- (Matrix::rankMatrix(WW) < nq)
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

#' Fit an MPST Model Using Distributed Learning
#'
#' @description Internal function to fit a multivariate penalized spline over 
#' triangulation (MPST) using distributed learning.
#'
#' @importFrom Matrix Matrix rankMatrix
#'
#' @param ind.Tr Index of the sub-triangulation.
#' @param Y The response variable vector.
#' @param Z A matrix of observation coordinates.
#' @param V A matrix of vertices of the triangulation.
#' @param Tr A triangulation matrix.
#' @param d The degree of piecewise polynomials.
#' @param r The smoothness parameter.
#' @param lambda The tuning parameter.
#' @param nl Number of lambda values to consider.
#' @param load.all Preloaded data for distributed processing.
#' @return A list containing distributed model fit components.
#' @keywords internal                          
fit.mpst.d <- function(ind.Tr, Y, Z, V, Tr, d = NULL, r = 1, lambda, nl, load.all) {
  # ns = parallel::detectCores()
  # TV = as.matrix(tdata(V, Tr)$TV)
  # inVT.list = inVT(V, Tr, Z)
  # load.all = worker.load(V = V, Tr = Tr, TV = TV, inVT.list = inVT.list, 
  #                        Y = Y, Z = Z, d = d, nl = nl, ns = ns)
  
  # 1. Preparation: extract data from subregions
  load.list = load.all[[ind.Tr]]
  
  Vs = load.list$Vs
  Trs = load.list$Trs
  inds = load.list$inds
  Ys = Y[inds]
  Zs = Z[inds, ]
  n = length(Ys)
  
  # 2. basis generation, smoothness conditions, and penalty function
  nd = ncol(Tr)
  if (nd == 3) {
    Bs <- as.matrix(basis2D.d(Vs, Trs, d, r, Zs)$B)
  } else if (nd == 4) {
    Bs <- as.matrix(basis3D.d(Vs, Trs, d, r, Zs)$B)
  }
  
  if (d < 1 | r < 0 | nrow(Trs) <= 1) {
    warning("The degree of Bernstein polynomials d has be greater than zero, the smoothness parameter r has be nonnegative!")
    H <- NA; Q2 <- NA;
  } else {
    if (nd == 3) {
      H <- as.matrix(smoothness2D(Vs, Trs, d, r))
    } else if (nd == 4) {
      H <- as.matrix(smoothness3D(Vs, Trs, d, r))
    }
    Q2 <- qrH(H)
  }
  
  if (d < 1) {
    warning("The degree of Bernstein polynomials d has be greater than zero.")
    K <- NA
  } else {
    if (nd == 3) {
      K <- as.matrix(energy(Vs, Trs, d))
    } else if (nd == 4) {
      K <- energy3D(Vs, Trs, d)
    }
  }
  
  # 2. Sub-model fitting
  ns <- length(Ys)
  nq.all <- ncol(Q2)
  
  W <- as.matrix(Bs %*% Q2)
  WW <- crossprod(W, W)
  rhs <- crossprod(W, Ys)
  D <- crossprod(t(crossprod(Q2, as.matrix(K))), Q2)
  D <- as.matrix(D)
  
  flag <- (Matrix::rankMatrix(WW) < nq.all)
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
    res <- Ys - beta
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
  res = Ys - beta
  mse = mean((Ys - Yhat)^2, na.rm = TRUE)
  
  nd = ncol(Tr)
  if (nd == 3) {
    nq = (d + 2) * (d + 1) / 2
  } else if (nd == 4) {
    nq = (d + 3) * (d + 2) * (d + 1) / 2 / 3
  }
  j = prodlim::row.match(load.list$t0, Trs)
  gamma.star = gamma[((j - 1) * nq + 1):(j * nq)]
  
  inVTs = inVT(Vs, Trs, Zs)
  ind.T1 = inVTs$ind.T
  B.star = Bs[which(ind.T1 == j), ((j - 1) * nq + 1):(j * nq)]

  mfit <- list(gamma.hat = gamma, 
               gamma.star = gamma.star,
               Bs = Bs,
               B.star = B.star,
               lambdac = lambdac,
               sse = sse,
               mse = mse,
               gcv = gcv,
               bic = bic, 
               edf = df,
               Vs = Vs,
               Trs = Trs)
  
  return(mfit)
}
