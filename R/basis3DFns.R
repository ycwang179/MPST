#' @importFrom Matrix sparseMatrix det
#' @importFrom pracma eps
#'
HQblkdiag3D <- function (A, cnt) { # Function HQblkdiag() should be the same for 2D & 3D
  n <- dim(A)[1]
  m <- dim(A)[2]
  k <- length(cnt) - 1
  D <- matrix(0, nrow = n, ncol = (k * m))
  for (i in 1:k) {
    D[(cnt[i] + 1):cnt[i + 1], ((i - 1) * m + 1):(i * m)] = A[(cnt[i] + 1):cnt[i + 1], ];
  }
  return(D)
}
# Example:
# A <- matrix(c(1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12), nrow = 4, ncol = 3)
# cnt <- matrix(c(1, 2, 3, 4), ncol = 1)
# HQblkdiag(A, cnt)

HQbary3D <- function (V, Th, X, Y, Z) { # Function HQbary() need to consider 3D locations instead of 2D
  k <- dim(Th)[1]
  A <- t(Th);
  B <- V[A, ];
  One <- matrix(1, nrow = dim(B)[1], ncol = 1)
  C <- cbind(One, B)
  cnt <- 4 * matrix(c(0:k), ncol = 1)
  C <- HQblkdiag(C, cnt)
  C <- t(C)

  X <- matrix(X, ncol = 1)
  Y <- matrix(Y, ncol = 1)
  Z <- matrix(Z, ncol = 1)
  One <- matrix(1, nrow = 1, ncol = dim(X)[1])
  D <- rbind(One, t(X), t(Y), t(Z))
  D <- t(matrix(rep(t(D), k), nrow = dim(t(D))[1]))
  Lam <- solve(C, D)
  return(Lam)
}
# Example:
# run eg1_main first
# X <- matrix(c(1, 2, 3), nrow = 1)
# Y <- matrix(c(4, 5, 6), nrow = 1)
# Z <- matrix(c(7, 8, 9), nrow = 1)
# HQbary3D(V, Th, X, Y, Z)

# HQgetInd() might be unnecessary;

loop3Dold <- function(d){
  nq <- choose(d + 3, 3)
  B <- matrix(0, nq, 4)
  l <- -1
  k <- -1
  a <- 0
  for (m in d:0) {
    l <- l + 1
    for (i in m:0) {
      k <- k + 1
      for (j in 0:i) {
        a <- a + 1
        B[a, ] <- c(i - j, j, k, l)
      }
    }
    k <- -1
  }
  return(B)
}

loop3D <- function (d) {
  tmp = expand.grid(d:0, 0:d, 0:d, 0:d)
  B = as.matrix(tmp[rowSums(tmp) == d, ])
  return(B)
}