#' Smoothness Constraints Matrix for Trivariate Spline over Tetrahedral Partition
#'
#' @importFrom Matrix sparseMatrix
#' @param V The \code{nV} by three matrix of vertices of a tetrahedron, where \code{nV} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The tetrahedral partition matrix of dimension \code{nT} by four, where \code{nT} is the number of tetrahedrons in the partition. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 9, and usually \code{d} is greater than one.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @return The smoothness constraints matrix.
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia.
#'
#' @examples
#' # example 1
#' d <- 4; r <- 1;
#' V <- matrix(rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 0), c(1, 1, 0),
#' c(0, 0, 1), c(0, 1, 1), c(1, 0, 1), c(1, 1, 1)), ncol = 3)
#' Tr <- matrix(rbind(c(5, 1, 2, 3), c(6, 5, 2, 3), c(6, 7, 5, 3),
#' c(6, 4, 7, 3), c(6, 2, 4, 3), c(6, 8, 7, 4)), ncol = 4)
#' H <- smoothness3D(V, Tr, d, r)
#' @export
#'
smoothness3D <- function (V, Tr, d, r) {
  Tr <- t(apply(t(Tr), 2, sort))
  l <- el3D(d, r)
  # t <- dtri(d)
  t = choose(d + 3, 3)
  thdata.list <- thdata(V, Tr)
  # Edge <- thdata.list$Edge
  Face <- thdata.list$Face
  TF <- thdata.list$TF.ind
  # TV <- thdata.list$TV
  # TE <- thdata.list$TE
  # FV <- thdata.list$FV
  # FE <- thdata.list$FE
  # EV <- thdata.list$EV
  # BF <- thdata.list$BF
  # BE <- thdata.list$BE
  # BV <- thdata.list$BV
  nT <- nrow(Tr)
  tf <- matrix(nrow = (max(TF[, 2])), ncol = 3)
  for (i in 1:(max(TF[, 2]))) {
    tf[i, 1] <- 1
    tf[i, 2] <- i
    tf[i, 3] <- length(which(TF[, 2] == i))
  }
  tf <- matrix(tf[, 3], nrow = 1)
  nF <- nrow(Face)
  tff <- tf - matrix(1, nrow = 1, ncol = nF)
  CF <- matrix(which(tff != 0), nrow = 1)
  ncf <- length(CF)
  SM.0 <- matrix(0, nrow = (ncf * l), ncol = (nT * t))
  SM <- matrix(ncol = 3)
  V <- t(V)
  TF.full <- sparseMatrix(i = TF[, 1], j = TF[, 2], x = TF[, 3])
  TF.full <- as.matrix(TF.full)
  for (i in 1:ncf) {
    Dum <- which(TF.full[, CF[1, i]] != 0, arr.ind = TRUE)
    t1 <- Dum[1]; t2 <- Dum[2]
    A <- Tr[t1, ]; B <- Tr[t2, ]
    A <- matrix(sort(A), nrow = 1); B <- matrix(sort(B), nrow = 1);
    T1 <- cbind(V[, A[1, 1]], V[, A[1, 2]], V[, A[1, 3]], V[, A[1, 4]])
    T2 <- cbind(V[, B[1, 1]], V[, B[1, 2]], V[, B[1, 3]], V[, B[1, 4]])
    H <- Lsmooth(T1, T2, d, r);
    SM.0[((i - 1) * l + 1):(i * l), ((t1 - 1) * t + 1):(t1 * t)] <- H[, 1:t]
    SM.0[((i - 1) * l + 1):(i * l), ((t2 - 1) * t + 1):(t2 * t)] <- H[, (t + 1):(2 * t)]
    }
  SM.1 <- which(SM.0 != 0, arr.ind = T)
  SM.Value <- matrix(nrow = dim(SM.1)[1], ncol = 1)
  for (j in 1:dim(SM.1)[1]) {
    SM.Value[j, 1] <- SM.0[SM.1[j, 1], SM.1[j, 2]]
  }
  SM <- cbind(SM.1, SM.Value)
  colnames(SM) <- c("i", "j", "x")
  
  # SM.full <- matrix(0, ncf * l, nT * t)
  # SM.initial <- as.matrix(sparseMatrix(i = SM[, 1], j = SM[, 2], x = SM[, 3]))
  # SM.full[1:dim(SM.initial)[1], 1:dim(SM.initial)[2]] <- SM.initial 
  # SM.sparse <- Matrix(SM.full, sparse = TRUE)
  
  SM.sparse <- Matrix(0, ncf * l, nT * t, sparse = TRUE)
  SM.initial <- sparseMatrix(i = SM[, 1], j = SM[, 2], x = SM[, 3])
  SM.sparse[1:dim(SM.initial)[1], 1:dim(SM.initial)[2]] <- SM.initial 
  
  return(SM.sparse)
}


######################################################################
# Other functions

# Lsmooth function
# This function produces a local smoothness matrix of a polynomial piece of degree d up
# to the order r on the tetrahedra which vertices are encoded in V and W.
Lsmooth <- function (V, W, d, r){
  comf.list.1 <- comf(V, W)
  comf.list.2 <- comf(W, V)
  CD <- cbind(V[, comf.list.1$C[1]], V[, comf.list.1$C[2]], V[, comf.list.1$C[3]], V[, comf.list.1$s], W[, comf.list.2$s])
  if ((length(comf.list.1$C) == 4) | (length(comf.list.1$C) < 3)) {
    H <- matrix(nrow = 0, ncol = 0)
  }
  X <- cbind(matrix(comf.list.1$C, nrow = 1), matrix(comf.list.1$D, nrow = 1))
  M <- matrix(c(rep(1, 12), rep(2, 4), rep(2, 8), rep(3, 8), rep(3, 4), rep(4, 12),
                rep(c(rep(1, 3), 2), 4), rep(c(2, 2, 3, 3), 4), rep(c(3, 4, 4, 4), 4)), nrow = 16, ncol = 6)
  k <- 0
  for (i in 1:16) {
    x <- matrix(M[i, ], 1)
    diff.X.x <- X-x
    if (all(diff.X.x == 0)) {
      k <- i
    }
  }
  if (k == 1) {
    H <- Scas1(CD, d, r)
  }else if (k == 2) {
    H <- Scas2(CD, d, r)
  }else if (k == 3) {
    H <- Scas3(CD, d, r)
  }else if (k == 4) {
    H <- Scas4(CD, d, r)
  }else if (k == 5) {
    H <- Scas5(CD, d, r)
  }else if (k == 6) {
    H <- Scas6(CD, d, r)
  }else if (k == 7) {
    H <- Scas7(CD, d, r)
  }else if (k == 8) {
    H <- Scas8(CD, d, r)
  }else if (k == 9) {
    H <- Scas9(CD, d, r)
  }else if (k == 10) {
    H <- Scas10(CD, d, r)
  }else if (k == 11) {
    H <- Scas11(CD, d, r)
  }else if (k == 12) {
    H <- Scas12(CD, d, r)
  }else if (k == 13) {
    H <- Scas13(CD, d, r)
  }else if (k == 14) {
    H <- Scas14(CD, d, r)
  }else if (k == 15) {
    H <- Scas15(CD, d, r)
  }else if (k == 16) {
    H <- Scas16(CD, d, r)
  }
  return(H)
}


# Scas1 function:
# This function stands for Smoothness case 1
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 1 the tetrahedra are <v1,v2,v3,v> and <v1,v2,v3,w>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas1 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose (d - m + 2, 2)
  }
  V <- W[, 1:4]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4]
    if (m <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(m)
      # g <- dtri(m)
      g = choose(m + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, j + beta, k + gamma, tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(m, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}

# Scas2 function:
# This function stands for Smoothness case 2
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 2 the tetrahedra are <v1,v2,v3,v> and <v1,v2,w,v3>. This suggests
# that an arrangement has been made in W
# H has size [l,n]
Scas2 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, 1:4]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (k <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(k)
      # g <- dtri(k)
      g = choose(k + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, j + beta, m + gamma, tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(k, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas3 function:
# This function stands for Smoothness case 3
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 3 the tetrahedra are <v1,v2,v3,v> and <v1,w,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas3 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, 1:4]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (j <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(j)
      # g <- dtri(j)
      g = choose(j + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, k + beta, m + gamma, tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(j, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas4 function:
# This function stands for Smoothness case 4
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 4 the tetrahedra are <v1,v2,v3,v> and <w,v1,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas4 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, 1:4]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (i <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(i)
      # g <- dtri(i)
      g = choose(i + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(j + alpha, k + beta, m + gamma, tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(i, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas5 function:
# This function stands for Smoothness case 5
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 5 the tetrahedra are <v1,v2,v,v3> and <v1,v2,v3,w>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas5 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 2, 4, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (m <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(m)
      # g <- dtri(m)
      g = choose(m + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, j + beta, gamma, k + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(m, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas6 function:
# This function stands for Smoothness case 6
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 6 the tetrahedra are <v1,v2,v,v3> and <v1,v2,w,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas6 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 2, 4, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (k <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(k)
      # g <- dtri(k)
      g = choose(k + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, j + beta, gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(k, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas7 function:
# This function stands for Smoothness case 7
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 7 the tetrahedra are <v1,v2,v,v3> and <v1,w,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas7 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 2, 4, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (j <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(j)
      # g <- dtri(j)
      g = choose(j + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, k + beta, gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(j, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas8 function:
# This function stands for Smoothness case 8
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 8 the tetrahedra are <v1,v2,v,v3> and <w,v1,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas8 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 2, 4, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (i <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(i)
      # g <- dtri(i)
      g = choose(i + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(j + alpha, k + beta, gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(i, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas9 function:
# This function stands for Smoothness case 9
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 9 the tetrahedra are <v1,v,v2,v3> and <v1,v2,v3,w>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas9 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 4, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4]
    if (m <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(m)
      # g <- dtri(m)
      g = choose(m + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, beta, j + gamma, k + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(m, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas10 function:
# This function stands for Smoothness case 10
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 9 the tetrahedra are <v1,v,v2,v3> and <v1,v2,w,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas10 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 4, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (k <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(k)
      # g <- dtri(k)
      g = choose(k + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, beta, j + gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(k, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas11 function:
# This function stands for Smoothness case 11
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 11 the tetrahedra are <v1,v,v2,v3> and <v1,w,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas11 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 4, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (j <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(j)
      # g <- dtri(j)
      g = choose(j + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(i + alpha, beta, k + gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(j, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas12 function:
# This function stands for Smoothness case 12
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 12 the tetrahedra are <v1,v,v2,v3> and <w,v1,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas12 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(1, 4, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (i <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(i)
      # g <- dtri(i)
      g = choose(i + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(j + alpha, beta, k + gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(i, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas13 function:
# This function stands for Smoothness case 13
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 12 the tetrahedra are <v,v1,v2,v3> and <v1,v2,v3,w>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas13 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(4, 1, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4]
    if (m <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(m)
      # g <- dtri(m)
      g = choose(m + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(alpha, i + beta, j + gamma, k + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(m, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas14 function:
# This function stands for Smoothness case 14
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 12 the tetrahedra are <v,v1,v2,v3> and <v1,v2,w,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas14 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(4, 1, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (k <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(k)
      # g <- dtri(k)
      g = choose(k + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(alpha, i + beta, j + gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(k, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas15 function:
# This function stands for Smoothness case 15
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 12 the tetrahedra are <v,v1,v2,v3> and <v1,w,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas15 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(4, 1, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i <- D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s,4];
    if (j <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(j)
      # g <- dtri(j)
      g = choose(j + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(alpha, i + beta, k + gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(j, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}


# Scas16 function:
# This function stands for Smoothness case 16
# W <- [v1,v2,v3,v,w] and we are interested in finding the smoothness
# conditions between <v1,v2,v3,v> and <v1,v2,v3,w>.
# in case 12 the tetrahedra are <v,v1,v2,v3> and <w,v1,v2,v3>. This assumes
# that some kind of arrangement has been made in W
# H has size [l,n]
Scas16 <- function (W, d, r){
  # f <- dtri(d)
  f = choose(d + 3, 3)
  n <- 2 * f
  l <- 0
  for (m in 0:r) {
    l <- l + choose(d - m + 2, 2)
  }
  V <- W[, c(4, 1, 2, 3)]
  Br <- bary3D(matrix(W[, 5], dim(W)[1], 1), V)
  H <- matrix(0, l, n)
  D <- loop3D(d)
  p <- 0
  for (s in 1:f){
    i=D[s, 1]; j <- D[s, 2]; k <- D[s, 3]; m <- D[s, 4];
    if (i <= r) {
      p <- p + 1
      q <- reperer(d, matrix(c(i, j, k, m), nrow = 1, ncol = 4))
      H[p, f + q] <- -1
      B <- loop3D(i)
      # g <- dtri(i)
      g = choose(i + 3, 3)
      for (t in 1:g){
        alpha <- B[t, 1]; beta <- B[t, 2]; gamma <- B[t, 3]; tau <- B[t, 4];
        re <- reperer(d, matrix(c(alpha, j + beta, k + gamma, m + tau), nrow = 1, ncol = 4))
        H[p, re] <- Bsn(i, matrix(c(alpha, beta, gamma, tau), nrow = 1), t(Br$Br))
      }
    }
  }
  return(H)
}
