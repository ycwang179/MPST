# Other Functions for Energy

# # 1. indices3d function: (# I think this is the same function as loop.)
# indices3d <- function (d) {
#   nq = choose(d + 3, 3)
#   B <- matrix(0, nrow = nq, ncol = 4)
#   l <- -1
#   k <- -1
#   a <- 0
#   for (m in d:0) {
#     l <- l + 1
#     for (i in m:0) {
#       k <- k + 1;
#       for (j in 0:i) {
#         a <- a + 1
#         B[a, ] <- c(i-j, j, k, l)
#       }
#     }
#     k <- -1
#   }
#   return(B)
# }
# # Example: indices3d(1)

# 2. chooseEngery function:
chooseEngery <- function (p, q) {
  if (is.matrix(p) & is.matrix(q)) {
    c <- matrix(0, nrow = dim(p)[1], ncol = dim(p)[2])
    aux <- matrix(0, nrow = dim(p)[1], ncol = dim(p)[2])
    for (i in 1:dim(p)[1]) {
      for (j in 1:dim(p)[2]) {
        c[i, j] <- factorial(q[i, j]) / (factorial(p[i, j]))
        aux[i, j] <- factorial(q[i, j] - p[i, j])
        if (aux[i, j] == 0) {
          c[i, j] <- 0
          aux[i, j] <- 1
        }
        c[i, j] <- c[i, j] / aux[i, j]
      }
    }
  } else if (is.numeric(p) & is.matrix(q)) {
    p <- matrix(p, nrow = dim(q)[1], ncol = dim(q)[2])
    c <- matrix(0, nrow = dim(p)[1], ncol = dim(p)[2])
    aux <- matrix(0, nrow = dim(p)[1], ncol = dim(p)[2])
    for (i in 1:dim(p)[1]) {
      for (j in 1:dim(p)[2]) {
        c[i, j] <- factorial(q[i, j]) / (factorial(p[i, j]))
        aux[i, j] <- factorial(q[i, j] - p[i, j])
        if (aux[i, j] == 0) {
          c[i, j] <- 0
          aux[i, j] <- 1
        }
        c[i, j] <- c[i, j] / aux[i, j]
      }
    }
  } else if (is.numeric(q) & is.matrix(p)) {
    q <- matrix(q, nrow = dim(p)[1], ncol = dim(p)[2])
    c <- matrix(0, nrow = dim(p)[1], ncol = dim(p)[2])
    aux <- matrix(0, nrow = dim(p)[1], ncol = dim(p)[2])
    for (i in 1:dim(p)[1]) {
      for (j in 1:dim(p)[2]) {
        c[i, j] <- factorial(q[i, j]) / (factorial(p[i, j]))
        aux[i, j] <- factorial(q[i, j] - p[i, j])
        if (aux[i, j] == 0) {
          c[i, j] <- 0
          aux[i, j] <- 1
        }
        c[i, j] <- c[i, j] / aux[i, j]
      }
    }
  } else if (is.numeric(q) & is.numeric(p)) {
    c <- factorial(q) / (factorial(p))
    aux <- factorial(q - p)
    I <- which(aux == 0)
    c[I] <- 0
    aux[I] <- 1
    c <- c / aux
  }
  return(c)
}
# Example:
# chooseEngery(3, 10)
# p <- matrix(c(1, 2, 1, 2), 2, 2)
# q <- matrix(c(6, 2, 8, 10), 2, 2)
# p <- 2
# q <- 3
# chooseEngery(p, q)


# 3. build function:
build3D <- function (d) {
  # B <- indices3d(d)
  B <- loop3D(d)
  # f <- dtri(d)
  f = choose(d + 3, 3)
  unit <- 1:f
  repunit1 <- t(matrix(rep(unit, f), nrow = f, ncol = f))
  R <- matrix(repunit1, ncol = 1)
  BR1 <- NULL
  for (i in 1:dim(B)[1]) {
    BR0 <- t(matrix(rep(B[i, ], f), nrow = 4, ncol = f))
    BR1 <- rbind(BR1, BR0)
  }
  BR2 <- NULL
  for (i in 1:f) {
    BR2 <- rbind(BR2, B)
  }
  BR <- BR1 + BR2
  G <- chooseEngery(BR1, BR)
  G <- matrix(apply(t(G), 2, prod), nrow = 1)
  G <- G / (6 * chooseEngery(d, 2 * d) * chooseEngery(3, 2 * d + 3))
  G <- matrix(G, f, f)
  return(G)
}
# Example: build3D(2)
