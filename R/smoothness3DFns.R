# Other Functions for Smoothness

# 1. Bsn function:
Bsn <- function(d, A, B){
  i <- A[1]; j <- A[2]; k <- A[3]; l <- A[4];
  m <- B[1]; n <- B[2]; p <- B[3]; q <- B[4];
  f <- factorial(d)/(factorial(i) * factorial(j) * factorial(k) * factorial(l)) * m^i * n^j * p^k * q^l
  return(f)
}
# Example:
# d <- 4
# A <- matrix(c(0, 0, 0, 0), nrow = 1)
# B <- matrix(c(1, 1, 1, 1), nrow = 1)
# Bsn(d, A, B)

# 2. comf function:
# This function checks if the tetrahedra which vertices are encoded in V
# and W share a common face and returns in A and B the indices of the
# vectors in the common face.
# If length of A or B is greater than 3, the tetrahedra are the same;
# If it's 3 then they have a common face;
# If it's 2 then they have a common edge;
# If it's 1, a common vertex;
# It it's 0, nothing in common.
# In V and W the columns give the coordinates of the vertices.
# s is the index of the vector of V not in W
comf <- function(V, W){
  A <- rep(0, 4);
  B <- rep(0, 4);
  # The sizes of A and B may change towards the end.
  for(i in 1:4){
    Z <- cbind(W[, 1] - V[, i], W[, 2] - V[, i], W[, 3] - V[, i], W[, 4] - V[, i])
    z <- zerov(Z)
    if(length(z) == 0){
      A[i] <- 0
      s = i
    }else{
      A[i] <- i
      B[i] <- z
    }
  }

  C <- A[which(A != 0)]
  D <- B[which(B != 0)]
  comf.list <- list(C = C, D = D, s = s)
  return(comf.list)
}
# Example:
#V <- matrix(c(2, 1, 2, 3, 3, 3, 4, 3, 4, 5, 3, 4), nrow = 3, ncol = 4)
#W <- matrix(c(1, 2, 3, 2, 5, 4, 2, 5, 5, 4, 3, 6), nrow = 3, ncol = 4)
#V <- matrix(c(2, 1, 2, 4, 3, 4, 4, 5, 4, 3, 4, 3, 3, 3, 4, 3, 4, 5, 3, 4), nrow = 5, ncol = 4)
#W <- matrix(c(2, 1, 2, 4, 3, 4, 4, 5, 4, 3, 5, 3, 3, 3, 4, 3, 4, 5, 3, 4), nrow = 5, ncol = 4)
#V <- matrix(c(1, 0, 0.5, 0.5, 0,1, 0.5, 0.5, 1, 1, 0, 1), nrow = 3, ncol = 4)
#W <- matrix(c(0.5, 0.5, 0.5, 1, 0, 0.5, 0.5, 0, 1, 0.5, 0.5, 1), nrow = 3, ncol = 4)
#comf(V, W)

# 3. del_reps function:
# This function deletes all repetitions from the list lold.
del_reps <- function(lold){
  if (!!length(lold)){# (length(lold) != 0)
    if (is.matrix(lold)){
      lold <- matrix(lold, nrow = max(dim(lold)), ncol = 1)
      lnew <- unique(lold)
    }else{
      lnew <- unique(lold)
      lold <- t(lold)
    }
  }else{
    lnew <- lold
  }
  return(lnew)
}
# Example:
#lold <- c(1, 2, 3, 1)
#lold <- matrix(c(1, 2, 2, 3, 2, 2, 3, 1, 6), nrow = 3, ncol = 3)
#del_reps(lold)

# 4. dtri function:
# This function calculate the dimensions of the trivariate polynomial space.
dtri <- function(d){
  f <- (d + 3) * (d + 2) * (d + 1) / 3 / 2
  return(f)
}
# Example: dtri(4)

# 5. el function:
el <- function(d, r){
  l <- 0
  for (m in 0:r){
    l <- l + choose(d - m + 2, 2)
  }
  return(l)
}
# Example: el(4, 1)


# 6. loop function:
# This function produces a matrix which rows give the indices in the order wanted of the
# B-coefficients on a tetrahedron
loop <- function(d){
  B <- matrix(0, dtri(d), 4)
  l <- -1
  k <- -1
  a <- 0
  for(m in d:0){
    l <- l + 1
    for(i in m:0){
      k <- k + 1
      for(j in 0:i){
        a <- a + 1
        B[a, ] <- c(i - j, j, k, l)
      }
    }
    k <- -1
  }
  return(B)
}
# Example:
# d <- 4
# loop(d)

# 7. reperer function
# This function finds the index corresponding to the domain point
# associated with A = c(i, j, k, l)
reperer <- function(d, A){
  i <- A[1]; j <- A[2]; k <- A[3]; l <- A[4];
  #Based on d, we introduce K
  L <- matrix(0, 1, d + 2)
  L[1] <- 1
  for(i in 2:(d + 1)){
    L[i] <- L[i-1] + i
  }
  L <- matrix(c(0, L[seq(d + 1, 1, by = -1)]), nrow = 1)
  K <- matrix(c(0, seq(d + 1 - l, 1, by = -1)), nrow = 1)
  f <- sum(L[1, 1:(l + 1)]) + sum(K[1, 1:(k + 1)]) + (j + 1)
  return(f)
}
# Example:
#A <- matrix(c(0, 0, 0, 0), nrow = 1)
#A <- matrix(c(1, 1, 1, 1), nrow = 1)
#A <- matrix(c(1, 1, 0, 1), nrow = 1)
#d <- 4
#reperer(d, A)

# 8. zerov function:
# This function finds the indices of the columns with zero elements in a matrix.
# The iteration is not necessary
zerov <- function(M){
  nr <- nrow(M); nc <- ncol(M);
  A <- which(apply(M == 0, 2, sum) == nr)
  return(A)
}
# Example:
# M <- matrix(c(1, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 0), nrow = 4, ncol = 3)
# zerov(M)
