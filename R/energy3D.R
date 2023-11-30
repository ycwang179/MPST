#' Energy Function for Trivariate Spline over Tetrahedral Partition
#'
#' @importFrom pracma cross kron
#' @importFrom Matrix sparseMatrix
#' @importFrom pryr mem_used
#' @param V The \code{nV} by three matrix of vertices of a tetrahedron, where \code{nV} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The tetrahedral partition matrix of dimension \code{nT} by four, where \code{nT} is the number of tetrahedrons in the partition. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 9, and usually \code{d} is greater than one.
#' \cr
#' @return The energy matrix.
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia.
#'
#' @examples
#' # example 1
#' d <- 4
#' V <- matrix(rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 0), c(1, 1, 0),
#' c(0, 0, 1), c(0, 1, 1), c(1, 0, 1), c(1, 1, 1)), ncol = 3)
#' Tr <- matrix(rbind(c(5, 1, 2, 3), c(6, 5, 2, 3), c(6, 7, 5, 3),
#' c(6, 4, 7, 3), c(6, 2, 4, 3), c(6, 8, 7, 4)), ncol = 4)
#' P <- energyM3D(V, Tr, d)
#' @export
#'
energy3D <- function (V, Tr, d) {
  n <- dim(Tr)[1]
  m <- choose((d + 3), 3)
  Vl <- V[t(Tr), ]
  idx4 <- matrix(c(1:n), ncol = 1) * 4
  idx3 <- idx4 - 1
  idx2 <- idx3 - 1
  idx1 <- idx2 - 1
  E1 <- Vl[idx1, ] - Vl[idx4, ]
  E2 <- Vl[idx2, ] - Vl[idx4, ]
  E3 <- Vl[idx3, ] - Vl[idx4, ]
  Volv <- matrix(abs(apply(E1 * cross(E2, E3), 1, sum)), ncol = 1)
  ridx1 <- t(matrix(rep(1:4, 3), 4, 3))
  cidx1 <- matrix(rep(1:3, 4), 3, 4)
  rvidx <- matrix(rep(matrix(ridx1, ncol = 1), n), ncol = 1)
  cvidx <- matrix(rep(matrix(rep(1:3, 4), 3, 4), n), ncol = 1)
  addrow <- 4 * matrix((floor((0:(12 * n - 1)) / 12)), ncol = 1)
  addcol <- 3 * matrix((floor((0:(12 * n - 1)) / 12)), ncol = 1)
  Bm1 <- matrix(0, nrow = dim(addrow)[1], ncol = 3)
  Bm1[, 1] <- rvidx + addrow
  Bm1[, 2] <- cvidx + addcol
  Bm1[, 3] <- matrix(t(Vl), ncol = 1)
  Bm1 <- Bm1[rowSums(Bm1 == 0) == 0, ]
  Bm1 <- Bm1[order(Bm1[, 2]), ]
  Bm1.sparse.0 <- sparseMatrix(i = Bm1[, 1], j = Bm1[, 2], x = Bm1[, 3])
  Bm1.sparse.0 = as.matrix(Bm1.sparse.0)
  # Bm1.sparse.0 <- t_shallow(Bm1.sparse.0)
  Bm1.sparse.0 <- t(Bm1.sparse.0)
  Bm1.sparse.1 <- which(Bm1.sparse.0 != 0, arr.ind = T)
  Bm1.sparse <- matrix(nrow = dim(Bm1.sparse.1)[1], ncol = 3)
  Bm1.sparse[, 1:2] <- Bm1.sparse.1[, 1:2]
  suppressMessages(Bm1.sparse[, 3] <- Bm1.sparse.0[Bm1.sparse.0!=0])
  Bm1 <- Bm1.sparse
  sumbary <- matrix(c(1, 1, 1, 1), nrow = 1)
  Idn.1 <- which(diag(n) != 0, arr.ind = T)
  Idn <- cbind(Idn.1, matrix(1, nrow = dim(Idn.1)[1], ncol = 1))

  # for kron(Idn,sumbary)
  Idn.full <- as.matrix(sparseMatrix(i = Idn[, 1], j = Idn[, 2], x = Idn[, 3]))
  Idn.S.0 <- kron(Idn.full, sumbary)
  Idn.S.1 <- which(Idn.S.0 != 0, arr.ind = T)
  Idn.S.2 <- matrix(nrow = dim(Idn.S.1)[1], ncol = 3)
  Idn.S.2[, 1:2] <- Idn.S.1[, 1:2]
  Idn.S.2[, 3] <- Idn.S.0[Idn.S.0!=0]
  Idn.S.sparse <- sparseMatrix(i = Idn.S.2[, 1], j = Idn.S.2[, 2], x = Idn.S.2[, 3])
  # end kron(Idn,sumbary)
  
  Bm1.sparse <- sparseMatrix(i = Bm1.sparse[, 1], j = Bm1.sparse[, 2], x = Bm1.sparse[, 3])
  B.sparse.0 <- as.matrix(rbind(Bm1.sparse, Idn.S.sparse))
  B.sparse.1 <- which(B.sparse.0 != 0, arr.ind = T)
  B <- matrix(nrow = dim(B.sparse.1)[1], ncol = 3)
  B[, 1:2] <- B.sparse.1[, 1:2]
  B[, 3] <- B.sparse.0[B.sparse.0!=0]
  
  vxv <- matrix(rep(c(1, 0, 0), n), ncol = 1)
  vyv <- matrix(rep(c(0, 1, 0), n), ncol = 1)
  vzv <- matrix(rep(c(0, 0, 1), n), ncol = 1)
  sv <- matrix(1, nrow = n, ncol = 1)
  vz <- matrix(0, nrow = (3 * n), ncol = 1)
  B.full <- as.matrix(sparseMatrix(i = B[, 1], j = B[, 2], x = B[, 3]))
  
  bvdx <- solve(B.full, rbind(vxv, sv))
  bvdy <- solve(B.full, rbind(vyv, sv))
  bvdz <- solve(B.full, rbind(vzv, sv))
  av <- solve(B.full, rbind(vz, sv))
  thdx <- bvdx - av
  thdy <- bvdy - av
  thdz <- bvdz - av
  
  Mat1 <- build3D(d - 2)
  # B <- indices3d(d)
  B <- loop3D(d)
  # C <- indices3d(d - 1)
  C <- loop3D(d - 1)
  I <- matrix(C[, 1], ncol = 1)
  J <- matrix(C[, 2], ncol = 1)
  K <- matrix(C[, 3], ncol = 1)
  L <- matrix(C[, 4], ncol = 1)
  C.1a <- cbind(I + 1, J, K, L)
  C.2a <- cbind(I, J + 1, K, L)
  C.3a <- cbind(I, J, K + 1, L)
  C.4a <- cbind(I, J, K, L + 1)
  
  Index1a = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.1a))))
  Index2a = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.2a))))
  Index3a = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.3a))))
  Index4a = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.4a))))
  
  # Index1a <- NULL
  # Index2a <- NULL
  # Index3a <- NULL
  # Index4a <- NULL
  # 
  # for (i in 1:dim(B)[1]) { # for Index1a=find(ismember(B,[I+1,J,K,L],'rows'));
  #   for (j in 1:dim(C.1a)[1]) {
  #     if (identical(B[i, ], C.1a[j, ])) {
  #       Index1a <- rbind(Index1a, i)
  #       break
  #     }
  #   }
  #   for (j in 1:dim(C.2a)[1]) {
  #     if (identical(B[i, ], C.2a[j, ])) {
  #       Index2a <- rbind(Index2a, i)
  #       break
  #     }
  #   }
  #   for (j in 1:dim(C.3a)[1]) {
  #     if (identical(B[i, ], C.3a[j, ])) {
  #       Index3a <- rbind(Index3a, i)
  #       break
  #     }
  #   }
  #   for (j in 1:dim(C.4a)[1]) {
  #     if (identical(B[i,], C.4a[j, ])) {
  #       Index4a <- rbind(Index4a, i)
  #       break
  #     }
  #   }
  # }
  # Index1a <- matrix(Index1a, ncol = 1)
  # Index2a <- matrix(Index2a, ncol = 1)
  # Index3a <- matrix(Index3a, ncol = 1)
  # Index4a <- matrix(Index4a, ncol = 1)
  
  Bin.1 <- which(diag(m) != 0, arr.ind = T)
  Bin <- cbind(Bin.1, matrix(1, nrow = dim(Bin.1)[1], ncol = 1))
  
  Bin.full <- as.matrix(sparseMatrix(i = Bin[, 1], j = Bin[, 2], x = Bin[, 3]))
  S1.0 <- Bin.full[Index1a, ]
  S1.1 <- which(S1.0 != 0, arr.ind = T)
  S1 <- cbind(S1.1, matrix(1, nrow = dim(S1.1)[1], ncol = 1))
  S2.0 <- Bin.full[Index2a, ]
  S2.1 <- which(S2.0 != 0, arr.ind = T)
  S2 <- cbind(S2.1, matrix(1, nrow = dim(S2.1)[1], ncol = 1))
  S3.0 <- Bin.full[Index3a, ]
  S3.1 <- which(S3.0 != 0, arr.ind = T)
  S3 <- cbind(S3.1, matrix(1, nrow = dim(S3.1)[1], ncol = 1))
  S4.0 <- Bin.full[Index4a, ]
  S4.1 <- which(S4.0 != 0, arr.ind = T)
  S4 <- cbind(S4.1, matrix(1, nrow = dim(S4.1)[1], ncol = 1))
  

  # for kron(thdx(idx1),S1)
  thdx.idx1 <- matrix(thdx[idx1], ncol = 1)
  S1.full <- as.matrix(sparseMatrix(i = S1[, 1], j = S1[, 2], x = S1[, 3]))
  thdx.idx1.S1.0 <- kron(thdx.idx1, S1.full)
  thdx.idx1.S1.1 <- which(thdx.idx1.S1.0 != 0, arr.ind = T)
  thdx.idx1.S1.2 <- matrix(nrow = dim(thdx.idx1.S1.1)[1], ncol = 3)
  thdx.idx1.S1.2[, 1:2] <- thdx.idx1.S1.1[, 1:2]
  thdx.idx1.S1.2[, 3] <- thdx.idx1.S1.0[thdx.idx1.S1.0!=0]
  if (nrow(thdx.idx1.S1.2)==0) {
    thdx.idx1.S1.full <- thdx.idx1.S1.0
  } else {
    thdx.idx1.S1.sparse <- sparseMatrix(i = thdx.idx1.S1.2[, 1], j = thdx.idx1.S1.2[, 2], x = thdx.idx1.S1.2[, 3])
    thdx.idx1.S1.full <- as.matrix(thdx.idx1.S1.sparse)
  }
  # end kron(thdx(idx1),S1)
  
  # for kron(thdx(idx2),S2)
  thdx.idx2 <- matrix(thdx[idx2], ncol = 1)
  S2.full <- as.matrix(sparseMatrix(i = S2[, 1], j = S2[, 2], x = S2[, 3]))
  thdx.idx2.S2.0 <- kron(thdx.idx2, S2.full)
  thdx.idx2.S2.1 <- which(thdx.idx2.S2.0 != 0, arr.ind = T)
  thdx.idx2.S2.2 <- matrix(nrow = dim(thdx.idx2.S2.1)[1], ncol = 3)
  thdx.idx2.S2.2[, 1:2] <- thdx.idx2.S2.1[, 1:2]
  thdx.idx2.S2.2[, 3] <- thdx.idx2.S2.0[thdx.idx2.S2.0!=0]
  if (nrow(thdx.idx2.S2.2)==0) {
    thdx.idx2.S2.full <- thdx.idx2.S2.0
  } else {
    thdx.idx2.S2.sparse <- sparseMatrix(i = thdx.idx2.S2.2[, 1], j = thdx.idx2.S2.2[, 2], x = thdx.idx2.S2.2[, 3])
    thdx.idx2.S2.full <- as.matrix(thdx.idx2.S2.sparse) 
  }
  # end kron(thdx(idx2),S2)
  
  # for kron(thdx(idx3),S3)
  thdx.idx3 <- matrix(thdx[idx3], ncol = 1)
  S3.full <- as.matrix(sparseMatrix(i = S3[, 1], j = S3[, 2], x = S3[, 3]))
  thdx.idx3.S3.0 <- kron(thdx.idx3, S3.full)
  thdx.idx3.S3.1 <- which(thdx.idx3.S3.0 != 0, arr.ind = T)
  thdx.idx3.S3.2 <- matrix(nrow = dim(thdx.idx3.S3.1)[1], ncol = 3)
  thdx.idx3.S3.2[, 1:2] <- thdx.idx3.S3.1[, 1:2]
  thdx.idx3.S3.2[, 3] <- thdx.idx3.S3.0[thdx.idx3.S3.0!=0]
  if (nrow(thdx.idx3.S3.2)==0) {
    thdx.idx3.S3.full <- thdx.idx3.S3.0
  } else {
    thdx.idx3.S3.sparse <- sparseMatrix(i = thdx.idx3.S3.2[, 1], j = thdx.idx3.S3.2[, 2], x = thdx.idx3.S3.2[, 3])
    thdx.idx3.S3.full <- as.matrix(thdx.idx3.S3.sparse)  
  }
  # end kron(thdx(idx3),S3)
  
  # for kron(thdx(idx4),S4)
  thdx.idx4 <- matrix(thdx[idx4], ncol = 1)
  S4.full <- as.matrix(sparseMatrix(i = S4[, 1], j = S4[, 2], x = S4[, 3]))
  thdx.idx4.S4.0 <- kron(thdx.idx4, S4.full)
  thdx.idx4.S4.1 <- which(thdx.idx4.S4.0 != 0, arr.ind = T)
  thdx.idx4.S4.2 <- matrix(nrow = dim(thdx.idx4.S4.1)[1], ncol = 3)
  thdx.idx4.S4.2[, 1:2] <- thdx.idx4.S4.1[, 1:2]
  thdx.idx4.S4.2[, 3] <- thdx.idx4.S4.0[thdx.idx4.S4.0!=0]
  if (nrow(thdx.idx4.S4.2)==0) {
    thdx.idx4.S4.full <- thdx.idx4.S4.0
  } else {
    thdx.idx4.S4.sparse <- sparseMatrix(i = thdx.idx4.S4.2[, 1], j = thdx.idx4.S4.2[, 2], x = thdx.idx4.S4.2[, 3])
    thdx.idx4.S4.full <- as.matrix(thdx.idx4.S4.sparse)
  }
  # end kron(thdx(idx4),S4)
  

  # for Dxgv
  thdx.idx1.S1.full.0 <- matrix(0, nrow = (dim(thdx.idx1)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdx.idx1.S1.full.0[1:dim(thdx.idx1.S1.full)[1], 1:dim(thdx.idx1.S1.full)[2]] <- thdx.idx1.S1.full
  
  thdx.idx2.S2.full.0 <- matrix(0, nrow = (dim(thdx.idx2)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdx.idx2.S2.full.0[1:dim(thdx.idx2.S2.full)[1], 1:dim(thdx.idx2.S2.full)[2]] <- thdx.idx2.S2.full
  
  thdx.idx3.S3.full.0 <- matrix(0, nrow = (dim(thdx.idx3)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdx.idx3.S3.full.0[1:dim(thdx.idx3.S3.full)[1], 1:dim(thdx.idx3.S3.full)[2]] <- thdx.idx3.S3.full
  
  thdx.idx4.S4.full.0 <- matrix(0, nrow = (dim(thdx.idx4)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdx.idx4.S4.full.0[1:dim(thdx.idx4.S4.full)[1], 1:dim(thdx.idx4.S4.full)[2]] <- thdx.idx4.S4.full
  
  thdx.idx.full.0 <- thdx.idx1.S1.full.0 + thdx.idx2.S2.full.0 + thdx.idx3.S3.full.0 + thdx.idx4.S4.full.0
  thdx.idx.full.1 <- which(thdx.idx.full.0 != 0, arr.ind = T)
  thdx.idx.full.2 <- matrix(nrow = dim(thdx.idx.full.1)[1], ncol = 3)
  thdx.idx.full.2[, 1:2] <- thdx.idx.full.1[, 1:2]
  thdx.idx.full.2[, 3] <- thdx.idx.full.0[thdx.idx.full.0!=0]
  thdx.idx.full.2[, 3] <- d * thdx.idx.full.2[, 3]
  Dxgv <- thdx.idx.full.2
  # end Dxgv
  
  # for kron(thdy(idx1),S1)
  thdy.idx1 <- matrix(thdy[idx1], ncol = 1)
  S1.full <- as.matrix(sparseMatrix(i = S1[, 1], j = S1[, 2], x = S1[, 3]))
  thdy.idx1.S1.0 <- kron(thdy.idx1, S1.full)
  thdy.idx1.S1.1 <- which(thdy.idx1.S1.0 != 0, arr.ind = T)
  thdy.idx1.S1.2 <- matrix(nrow = dim(thdy.idx1.S1.1)[1], ncol = 3)
  thdy.idx1.S1.2[, 1:2] <- thdy.idx1.S1.1[, 1:2]
  thdy.idx1.S1.2[, 3] <- thdy.idx1.S1.0[thdy.idx1.S1.0!=0]
  if (nrow(thdy.idx1.S1.2)==0) {
    thdy.idx1.S1.full <- thdy.idx1.S1.0
  } else {
    thdy.idx1.S1.sparse <- sparseMatrix(i = thdy.idx1.S1.2[, 1], j = thdy.idx1.S1.2[, 2], x = thdy.idx1.S1.2[, 3])
    thdy.idx1.S1.full <- as.matrix(thdy.idx1.S1.sparse) 
  }
  # end kron(thdy(idx1),S1)
  
  # for kron(thdy(idx2),S2)
  thdy.idx2 <- matrix(thdy[idx2], ncol = 1)
  S2.full <- as.matrix(sparseMatrix(i = S2[, 1], j = S2[, 2], x = S2[, 3]))
  thdy.idx2.S2.0 <- kron(thdy.idx2, S2.full)
  thdy.idx2.S2.1 <- which(thdy.idx2.S2.0 != 0, arr.ind = T)
  thdy.idx2.S2.2 <- matrix(nrow = dim(thdy.idx2.S2.1)[1], ncol = 3)
  thdy.idx2.S2.2[, 1:2] <- thdy.idx2.S2.1[, 1:2]
  thdy.idx2.S2.2[, 3] <- thdy.idx2.S2.0[thdy.idx2.S2.0!=0]
  if (nrow(thdy.idx2.S2.2)==0) {
    thdy.idx2.S2.full <- thdy.idx2.S2.0
  } else {
    thdy.idx2.S2.sparse <- sparseMatrix(i = thdy.idx2.S2.2[, 1], j = thdy.idx2.S2.2[, 2], x = thdy.idx2.S2.2[, 3])
    thdy.idx2.S2.full <- as.matrix(thdy.idx2.S2.sparse) 
  } 
  # end kron(thdy(idx2),S2)
  
  # for kron(thdy(idx3),S3)
  thdy.idx3 <- matrix(thdy[idx3], ncol = 1)
  S3.full <- as.matrix(sparseMatrix(i = S3[, 1], j = S3[, 2], x = S3[, 3]))
  thdy.idx3.S3.0 <- kron(thdy.idx3, S3.full)
  thdy.idx3.S3.1 <- which(thdy.idx3.S3.0 != 0, arr.ind = T)
  thdy.idx3.S3.2 <- matrix(nrow = dim(thdy.idx3.S3.1)[1], ncol = 3)
  thdy.idx3.S3.2[, 1:2] <- thdy.idx3.S3.1[, 1:2]
  thdy.idx3.S3.2[, 3] <- thdy.idx3.S3.0[thdy.idx3.S3.0!=0]
  if (nrow(thdy.idx3.S3.2)==0) {
    thdy.idx3.S3.full <- thdy.idx3.S3.0
  } else {
    thdy.idx3.S3.sparse <- sparseMatrix(i = thdy.idx3.S3.2[, 1], j = thdy.idx3.S3.2[, 2], x = thdy.idx3.S3.2[, 3])
    thdy.idx3.S3.full <- as.matrix(thdy.idx3.S3.sparse)
  } 
  # end kron(thdy(idx3),S3)
  
  # for kron(thdy(idx4),S4)
  thdy.idx4 <- matrix(thdy[idx4], ncol = 1)
  S4.full <- as.matrix(sparseMatrix(i = S4[, 1], j = S4[, 2], x = S4[, 3]))
  thdy.idx4.S4.0 <- kron(thdy.idx4, S4.full)
  thdy.idx4.S4.1 <- which(thdy.idx4.S4.0 != 0, arr.ind = T)
  thdy.idx4.S4.2 <- matrix(nrow = dim(thdy.idx4.S4.1)[1], ncol = 3)
  thdy.idx4.S4.2[, 1:2] <- thdy.idx4.S4.1[, 1:2]
  thdy.idx4.S4.2[, 3] <- thdy.idx4.S4.0[thdy.idx4.S4.0!=0]
  if (nrow(thdy.idx4.S4.2)==0) {
    thdy.idx4.S4.full <- thdy.idx4.S4.0
  } else {
    thdy.idx4.S4.sparse <- sparseMatrix(i = thdy.idx4.S4.2[, 1], j = thdy.idx4.S4.2[, 2], x = thdy.idx4.S4.2[, 3])
    thdy.idx4.S4.full <- as.matrix(thdy.idx4.S4.sparse) 
  } 
  # end kron(thdy(idx4),S4)
  

  # for Dygv
  thdy.idx1.S1.full.0 <- matrix(0, nrow = (dim(thdy.idx1)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdy.idx1.S1.full.0[1:dim(thdy.idx1.S1.full)[1], 1:dim(thdy.idx1.S1.full)[2]] <- thdy.idx1.S1.full
  
  thdy.idx2.S2.full.0 <- matrix(0, nrow = (dim(thdy.idx2)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdy.idx2.S2.full.0[1:dim(thdy.idx2.S2.full)[1], 1:dim(thdy.idx2.S2.full)[2]] <- thdy.idx2.S2.full
  
  thdy.idx3.S3.full.0 <- matrix(0, nrow = (dim(thdy.idx3)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdy.idx3.S3.full.0[1:dim(thdy.idx3.S3.full)[1], 1:dim(thdy.idx3.S3.full)[2]] <- thdy.idx3.S3.full
  
  thdy.idx4.S4.full.0 <- matrix(0, nrow = (dim(thdy.idx4)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdy.idx4.S4.full.0[1:dim(thdy.idx4.S4.full)[1], 1:dim(thdy.idx4.S4.full)[2]] <- thdy.idx4.S4.full
  
  thdy.idx.full.0 <- thdy.idx1.S1.full.0 + thdy.idx2.S2.full.0 + thdy.idx3.S3.full.0 + thdy.idx4.S4.full.0
  thdy.idx.full.1 <- which(thdy.idx.full.0 != 0, arr.ind = T)
  thdy.idx.full.2 <- matrix(nrow = dim(thdy.idx.full.1)[1], ncol = 3)
  thdy.idx.full.2[, 1:2] <- thdy.idx.full.1[, 1:2]
  thdy.idx.full.2[, 3] <- thdy.idx.full.0[thdy.idx.full.0!=0]
  thdy.idx.full.2[, 3] <- d * thdy.idx.full.2[, 3]
  Dygv <- thdy.idx.full.2
  # end Dygv
  
  # for kron(thdz(idx1),S1)
  thdz.idx1 <- matrix(thdz[idx1], ncol = 1)
  S1.full <- as.matrix(sparseMatrix(i = S1[, 1], j = S1[, 2], x = S1[, 3]))
  thdz.idx1.S1.0 <- kron(thdz.idx1, S1.full)
  thdz.idx1.S1.1 <- which(thdz.idx1.S1.0 != 0, arr.ind = T)
  thdz.idx1.S1.2 <- matrix(nrow = dim(thdz.idx1.S1.1)[1], ncol = 3)
  thdz.idx1.S1.2[, 1:2] <- thdz.idx1.S1.1[, 1:2]
  thdz.idx1.S1.2[, 3] <- thdz.idx1.S1.0[thdz.idx1.S1.0!=0]
  if (nrow(thdz.idx1.S1.2)==0) {
    thdz.idx1.S1.full <- thdz.idx1.S1.0
  } else {
    thdz.idx1.S1.sparse <- sparseMatrix(i = thdz.idx1.S1.2[, 1], j = thdz.idx1.S1.2[, 2], x = thdz.idx1.S1.2[, 3])
    thdz.idx1.S1.full <- as.matrix(thdz.idx1.S1.sparse) 
  }
  # end kron(thdz(idx1),S1)
  
  # for kron(thdz(idx2),S2)
  thdz.idx2 <- matrix(thdz[idx2], ncol = 1)
  S2.full <- as.matrix(sparseMatrix(i = S2[, 1], j = S2[, 2], x = S2[, 3]))
  thdz.idx2.S2.0 <- kron(thdz.idx2, S2.full)
  thdz.idx2.S2.1 <- which(thdz.idx2.S2.0 != 0, arr.ind = T)
  thdz.idx2.S2.2 <- matrix(nrow = dim(thdz.idx2.S2.1)[1], ncol = 3)
  thdz.idx2.S2.2[, 1:2] <- thdz.idx2.S2.1[, 1:2]
  thdz.idx2.S2.2[, 3] <- thdz.idx2.S2.0[thdz.idx2.S2.0!=0]
  if (nrow(thdz.idx2.S2.2)==0) {
    thdz.idx2.S2.full <- thdz.idx2.S2.0
  } else {
    thdz.idx2.S2.sparse <- sparseMatrix(i = thdz.idx2.S2.2[, 1], j = thdz.idx2.S2.2[, 2], x = thdz.idx2.S2.2[, 3])
    thdz.idx2.S2.full <- as.matrix(thdz.idx2.S2.sparse) 
  } 
  # end kron(thdz(idx2),S2)
  
  # for kron(thdz(idx3),S3)
  thdz.idx3 <- matrix(thdz[idx3], ncol = 1)
  S3.full <- as.matrix(sparseMatrix(i = S3[, 1], j = S3[, 2], x = S3[, 3]))
  thdz.idx3.S3.0 <- kron(thdz.idx3, S3.full)
  thdz.idx3.S3.1 <- which(thdz.idx3.S3.0 != 0, arr.ind = T)
  thdz.idx3.S3.2 <- matrix(nrow = dim(thdz.idx3.S3.1)[1], ncol = 3)
  thdz.idx3.S3.2[, 1:2] <- thdz.idx3.S3.1[, 1:2]
  thdz.idx3.S3.2[, 3] <- thdz.idx3.S3.0[thdz.idx3.S3.0!=0]
  if (nrow(thdz.idx3.S3.2)==0) {
    thdz.idx3.S3.full <- thdz.idx3.S3.0
  } else {
    thdz.idx3.S3.sparse <- sparseMatrix(i = thdz.idx3.S3.2[, 1], j = thdz.idx3.S3.2[, 2], x = thdz.idx3.S3.2[, 3])
    thdz.idx3.S3.full <- as.matrix(thdz.idx3.S3.sparse)
  }
  # end kron(thdz(idx3),S3)
  
  # for kron(thdz(idx4),S4)
  thdz.idx4 <- matrix(thdz[idx4], ncol = 1)
  S4.full <- as.matrix(sparseMatrix(i = S4[, 1], j = S4[,2], x = S4[, 3]))
  thdz.idx4.S4.0 <- kron(thdz.idx4, S4.full)
  thdz.idx4.S4.1 <- which(thdz.idx4.S4.0 != 0, arr.ind = T)
  thdz.idx4.S4.2 <- matrix(nrow = dim(thdz.idx4.S4.1)[1], ncol = 3)
  thdz.idx4.S4.2[, 1:2] <- thdz.idx4.S4.1[, 1:2]
  thdz.idx4.S4.2[, 3] <- thdz.idx4.S4.0[thdz.idx4.S4.0!=0]
  if (nrow(thdz.idx4.S4.2)==0) {
    thdz.idx4.S4.full <- thdz.idx4.S4.0
  } else {
    thdz.idx4.S4.sparse <- sparseMatrix(i=thdz.idx4.S4.2[, 1], j = thdz.idx4.S4.2[, 2], x = thdz.idx4.S4.2[, 3])
    thdz.idx4.S4.full <- as.matrix(thdz.idx4.S4.sparse) 
  } 
  # end kron(thdz(idx4),S4)
  

  # for Dzgv
  thdz.idx1.S1.full.0 <- matrix(0, nrow = (dim(thdz.idx1)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdz.idx1.S1.full.0[1:dim(thdz.idx1.S1.full)[1], 1:dim(thdz.idx1.S1.full)[2]] <- thdz.idx1.S1.full
  
  thdz.idx2.S2.full.0 <- matrix(0, nrow = (dim(thdz.idx2)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdz.idx2.S2.full.0[1:dim(thdz.idx2.S2.full)[1], 1:dim(thdz.idx2.S2.full)[2]] <- thdz.idx2.S2.full
  
  thdz.idx3.S3.full.0 <- matrix(0, nrow = (dim(thdz.idx3)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdz.idx3.S3.full.0[1:dim(thdz.idx3.S3.full)[1], 1:dim(thdz.idx3.S3.full)[2]] <- thdz.idx3.S3.full
  
  thdz.idx4.S4.full.0 <- matrix(0, nrow = (dim(thdz.idx4)[1] * dim(C)[1]), ncol = dim(Bin)[1])
  thdz.idx4.S4.full.0[1:dim(thdz.idx4.S4.full)[1], 1:dim(thdz.idx4.S4.full)[2]] <- thdz.idx4.S4.full
  
  thdz.idx.full.0 <- thdz.idx1.S1.full.0 + thdz.idx2.S2.full.0 + thdz.idx3.S3.full.0 + thdz.idx4.S4.full.0
  thdz.idx.full.1 <- which(thdz.idx.full.0 != 0, arr.ind = T)
  thdz.idx.full.2 <- matrix(nrow = dim(thdz.idx.full.1)[1], ncol = 3)
  thdz.idx.full.2[, 1:2] <- thdz.idx.full.1[, 1:2]
  thdz.idx.full.2[, 3] <- thdz.idx.full.0[thdz.idx.full.0!=0]
  thdz.idx.full.2[, 3] <- d * thdz.idx.full.2[, 3]
  Dzgv <- thdz.idx.full.2
  # end Dzgv
  
  # B <- indices3d(d - 1)
  B <- loop3D(d - 1)
  # C <- indices3d(d - 2)
  C <- loop3D(d - 2)
  I <- matrix(C[, 1], ncol = 1)
  J <- matrix(C[, 2], ncol = 1)
  K <- matrix(C[, 3], ncol = 1)
  L <- matrix(C[, 4], ncol = 1)
  C.1 <- cbind(I + 1, J, K, L)
  C.2 <- cbind(I, J + 1, K, L)
  C.3 <- cbind(I, J, K + 1, L)
  C.4 <- cbind(I, J, K, L + 1)
  
  Index1 = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.1))))
  Index2 = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.2))))
  Index3 = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.3))))
  Index4 = which(!is.na(prodlim::row.match(data.frame(B), data.frame(C.4))))
  
  # Index1 <- NULL
  # Index2 <- NULL
  # Index3 <- NULL
  # Index4 <- NULL
  # for (i in 1:dim(B)[1]) { # for Index1=find(ismember(B,[I+1,J,K,L],'rows'));
  #   for (j in 1:dim(C.1)[1]) {
  #     if (identical(B[i, ], C.1[j, ])) {
  #       Index1 <- rbind(Index1, i)
  #       break
  #     }
  #   }
  #   for (j in 1:dim(C.2)[1]) {
  #     if (identical(B[i, ], C.2[j, ])) {
  #       Index2 <- rbind(Index2, i)
  #       break
  #     }
  #   }
  #   for (j in 1:dim(C.3)[1]) {
  #     if (identical(B[i, ], C.3[j, ])) {
  #       Index3 <- rbind(Index3, i)
  #       break
  #     }
  #   }
  #   for (j in 1:dim(C.4)[1]) {
  #     if (identical(B[i, ], C.4[j, ])) {
  #       Index4 <- rbind(Index4, i)
  #       break
  #     }
  #   }
  # }
  # Index1 <- matrix(Index1, ncol = 1)
  # Index2 <- matrix(Index2, ncol = 1)
  # Index3 <- matrix(Index3, ncol = 1)
  # Index4 <- matrix(Index4, ncol = 1)
  
  sidx <- t(matrix(rep((0:(n-1)) * length(Index1a), length(Index1)), nrow = n))
  Ir1 <- matrix(rep(Index1, n), ncol = 1) + matrix(sidx, ncol = 1)
  Ir2 <- matrix(rep(Index2, n), ncol = 1) + matrix(sidx, ncol = 1)
  Ir3 <- matrix(rep(Index3, n), ncol = 1) + matrix(sidx, ncol = 1)
  Ir4 <- matrix(rep(Index4, n), ncol = 1) + matrix(sidx, ncol = 1)
  
  thxv1 <- t(matrix(rep(thdx.idx1, length(Index1)), nrow = length(thdx.idx1)))
  thxv2 <- t(matrix(rep(thdx.idx2, length(Index1)), nrow = length(thdx.idx2)))
  thxv3 <- t(matrix(rep(thdx.idx3, length(Index1)), nrow = length(thdx.idx3)))
  thxv4 <- t(matrix(rep(thdx.idx4, length(Index1)), nrow = length(thdx.idx4)))
  
  thyv1 <- t(matrix(rep(thdy.idx1, length(Index1)), nrow = length(thdy.idx1)))
  thyv2 <- t(matrix(rep(thdy.idx2, length(Index1)), nrow = length(thdy.idx2)))
  thyv3 <- t(matrix(rep(thdy.idx3, length(Index1)), nrow = length(thdy.idx3)))
  thyv4 <- t(matrix(rep(thdy.idx4, length(Index1)), nrow = length(thdy.idx4)))
  
  thzv1 <- t(matrix(rep(thdz.idx1, length(Index1)), nrow = length(thdz.idx1)))
  thzv2 <- t(matrix(rep(thdz.idx2, length(Index1)), nrow = length(thdz.idx2)))
  thzv3 <- t(matrix(rep(thdz.idx3, length(Index1)), nrow = length(thdz.idx3)))
  thzv4 <- t(matrix(rep(thdz.idx4, length(Index1)), nrow = length(thdz.idx4)))
  
  # for dxxv
  Dxgv.full <- as.matrix(sparseMatrix(i = Dxgv[, 1], j = Dxgv[, 2], x = Dxgv[, 3]))
  Dxgv.Ir1.thxv1 <- sweep(Dxgv.full[Ir1, ], 1, matrix(thxv1, ncol = 1), "*")
  Dxgv.Ir2.thxv2 <- sweep(Dxgv.full[Ir2, ], 1, matrix(thxv2, ncol = 1), "*")
  Dxgv.Ir3.thxv3 <- sweep(Dxgv.full[Ir3, ], 1, matrix(thxv3, ncol = 1), "*")
  Dxgv.Ir4.thxv4 <- sweep(Dxgv.full[Ir4, ], 1, matrix(thxv4, ncol = 1), "*")
  
  dxxv.0 <- (d - 1) * (Dxgv.Ir1.thxv1 + Dxgv.Ir2.thxv2 + Dxgv.Ir3.thxv3 + Dxgv.Ir4.thxv4)
  dxxv.1 <- which(dxxv.0 != 0, arr.ind = T)
  dxxv.2 <- matrix(nrow = dim(dxxv.1)[1], ncol = 3)
  dxxv.2[, 1:2] <- dxxv.1[, 1:2]
  dxxv.2[, 3] <- dxxv.0[dxxv.0!=0]
  dxxv <- dxxv.2
  # end dxxv
  
  # for dyyv
  Dygv.full <- as.matrix(sparseMatrix(i = Dygv[, 1], j = Dygv[, 2], x = Dygv[, 3]))
  Dygv.Ir1.thyv1 <- sweep(Dygv.full[Ir1, ], 1, matrix(thyv1, ncol = 1), "*")
  Dygv.Ir2.thyv2 <- sweep(Dygv.full[Ir2, ], 1, matrix(thyv2, ncol = 1), "*")
  Dygv.Ir3.thyv3 <- sweep(Dygv.full[Ir3, ], 1, matrix(thyv3, ncol = 1), "*")
  Dygv.Ir4.thyv4 <- sweep(Dygv.full[Ir4, ], 1, matrix(thyv4, ncol = 1), "*")
  
  dyyv.0 <- (d - 1) * (Dygv.Ir1.thyv1 + Dygv.Ir2.thyv2 + Dygv.Ir3.thyv3 + Dygv.Ir4.thyv4)
  dyyv.1 <- which(dyyv.0 != 0, arr.ind = T)
  dyyv.2 <- matrix(nrow = dim(dyyv.1)[1], ncol = 3)
  dyyv.2[, 1:2] <- dyyv.1[, 1:2]
  dyyv.2[, 3] <- dyyv.0[dyyv.0 != 0]
  dyyv <- dyyv.2
  # end dyyv
  
  # for dzzv
  Dzgv.full <- as.matrix(sparseMatrix(i = Dzgv[, 1], j = Dzgv[, 2], x = Dzgv[, 3]))
  Dzgv.Ir1.thzv1 <- sweep(Dzgv.full[Ir1, ], 1, matrix(thzv1, ncol = 1), "*")
  Dzgv.Ir2.thzv2 <- sweep(Dzgv.full[Ir2, ], 1, matrix(thzv2, ncol = 1), "*")
  Dzgv.Ir3.thzv3 <- sweep(Dzgv.full[Ir3, ], 1, matrix(thzv3, ncol = 1), "*")
  Dzgv.Ir4.thzv4 <- sweep(Dzgv.full[Ir4, ], 1, matrix(thzv4, ncol = 1), "*")
  
  dzzv.0 <- (d - 1) * (Dzgv.Ir1.thzv1 + Dzgv.Ir2.thzv2 + Dzgv.Ir3.thzv3 + Dzgv.Ir4.thzv4)
  dzzv.1 <- which(dzzv.0 != 0, arr.ind = T)
  dzzv.2 <- matrix(nrow = dim(dzzv.1)[1], ncol = 3)
  dzzv.2[, 1:2] <- dzzv.1[, 1:2]
  dzzv.2[, 3] <- dzzv.0[dzzv.0!=0]
  dzzv <- dzzv.2
  # end dzzv
  
  # for dxyv
  Dxgv.full <- as.matrix(sparseMatrix(i = Dxgv[, 1], j = Dxgv[, 2], x = Dxgv[, 3]))
  Dxgv.Ir1.thyv1 <- sweep(Dxgv.full[Ir1, ], 1, matrix(thyv1, ncol = 1), "*")
  Dxgv.Ir2.thyv2 <- sweep(Dxgv.full[Ir2, ], 1, matrix(thyv2, ncol = 1), "*")
  Dxgv.Ir3.thyv3 <- sweep(Dxgv.full[Ir3, ], 1, matrix(thyv3, ncol = 1), "*")
  Dxgv.Ir4.thyv4 <- sweep(Dxgv.full[Ir4, ], 1, matrix(thyv4, ncol = 1), "*")
  
  dxyv.0 <- (d - 1) * (Dxgv.Ir1.thyv1 + Dxgv.Ir2.thyv2 + Dxgv.Ir3.thyv3 + Dxgv.Ir4.thyv4)
  dxyv.1 <- which(dxyv.0 != 0, arr.ind = T)
  dxyv.2 <- matrix(nrow = dim(dxyv.1)[1], ncol = 3)
  dxyv.2[, 1:2] <- dxyv.1[, 1:2]
  dxyv.2[, 3] <- dxyv.0[dxyv.0!=0]
  dxyv <- dxyv.2
  # end dxyv
  
  # for dxzv
  Dxgv.full <- as.matrix(sparseMatrix(i = Dxgv[, 1], j = Dxgv[, 2], x = Dxgv[, 3]))
  Dxgv.Ir1.thzv1 <- sweep(Dxgv.full[Ir1, ], 1, matrix(thzv1, ncol = 1), "*")
  Dxgv.Ir2.thzv2 <- sweep(Dxgv.full[Ir2, ], 1, matrix(thzv2, ncol = 1), "*")
  Dxgv.Ir3.thzv3 <- sweep(Dxgv.full[Ir3, ], 1, matrix(thzv3, ncol = 1), "*")
  Dxgv.Ir4.thzv4 <- sweep(Dxgv.full[Ir4, ], 1, matrix(thzv4, ncol = 1), "*")
  
  dxzv.0 <- (d - 1) * (Dxgv.Ir1.thzv1 + Dxgv.Ir2.thzv2 + Dxgv.Ir3.thzv3 + Dxgv.Ir4.thzv4)
  dxzv.1 <- which(dxzv.0 != 0, arr.ind = T)
  dxzv.2 <- matrix(nrow = dim(dxzv.1)[1], ncol = 3)
  dxzv.2[, 1:2] <- dxzv.1[, 1:2]
  dxzv.2[, 3] <- dxzv.0[dxzv.0!=0]
  dxzv <- dxzv.2
  # end dxzv
  
  # for dyzv
  Dygv.full <- as.matrix(sparseMatrix(i = Dygv[, 1], j = Dygv[, 2], x = Dygv[, 3]))
  Dygv.Ir1.thzv1 <- sweep(Dygv.full[Ir1, ], 1, matrix(thzv1, ncol = 1), "*")
  Dygv.Ir2.thzv2 <- sweep(Dygv.full[Ir2, ], 1, matrix(thzv2, ncol = 1), "*")
  Dygv.Ir3.thzv3 <- sweep(Dygv.full[Ir3, ], 1, matrix(thzv3, ncol = 1), "*")
  Dygv.Ir4.thzv4 <- sweep(Dygv.full[Ir4, ], 1, matrix(thzv4, ncol = 1), "*")
  
  dyzv.0 <- (d - 1) * (Dygv.Ir1.thzv1 + Dygv.Ir2.thzv2 + Dygv.Ir3.thzv3 + Dygv.Ir4.thzv4)
  dyzv.1 <- which(dyzv.0 != 0, arr.ind = T)
  dyzv.2 <- matrix(nrow = dim(dyzv.1)[1], ncol = 3)
  dyzv.2[, 1:2] <- dyzv.1[, 1:2]
  dyzv.2[, 3] <- dyzv.0[dyzv.0!=0]
  dyzv <- dyzv.2
  # end dyzv
  
  cidx <- t(matrix(rep((1:m), (length(Index1)) * n), nrow = m))
  cs <- matrix(rep(t(matrix(rep((0:(n-1)) * m, length(Index1)[1]), nrow = n)), m), ncol = (n * m))
  dxxv.full <- as.matrix(sparseMatrix(i = dxxv[, 1], j = dxxv[, 2], x = dxxv[, 3]))
  ridx <- matrix(rep((1:dim(dxxv.full)[1]), dim(dxxv.full)[2]), ncol = 1)

  # for Dxx
  Dxx <- matrix(nrow = dim(ridx)[1], ncol = 3)
  Dxx[, 1] <- ridx[, 1]
  Dxx[, 2] <- matrix(cidx, ncol = 1) + matrix(cs, ncol = 1)
  Dxx[, 3] <- matrix(dxxv.full, ncol = 1)
  Dxx.0 <- sparseMatrix(i = Dxx[, 1], j = Dxx[, 2], x = Dxx[, 3])
  Dxx.0 = as.matrix(Dxx.0)
  Dxx.1 <- which(Dxx.0 != 0, arr.ind = T)
  Dxx <- matrix(nrow = dim(Dxx.1)[1], ncol = 3)
  Dxx[, 1:2] <- Dxx.1[, 1:2]
  suppressMessages(Dxx[, 3] <- Dxx.0[Dxx.0 != 0])
  # End Dxx
  
  # for Dyy
  dyyv.full <- as.matrix(sparseMatrix(i = dyyv[, 1], j = dyyv[, 2], x = dyyv[, 3]))
  Dyy <- matrix(nrow = dim(ridx)[1], ncol = 3)
  Dyy[, 1] <- ridx[, 1]
  Dyy[, 2] <- matrix(cidx, ncol = 1) + matrix(cs, ncol = 1)
  Dyy[, 3] <- matrix(dyyv.full, ncol = 1)
  Dyy.0 <- sparseMatrix(i = Dyy[, 1], j = Dyy[, 2], x = Dyy[, 3])
  Dyy.0 = as.matrix(Dyy.0)
  Dyy.1 <- which(Dyy.0 != 0, arr.ind = T)
  Dyy <- matrix(nrow = dim(Dyy.1)[1], ncol = 3)
  Dyy[, 1:2] <- Dyy.1[, 1:2]
  suppressMessages(Dyy[, 3] <- Dyy.0[Dyy.0!=0])
  # End Dyy
  
  # for Dzz
  dzzv.full <- as.matrix(sparseMatrix(i = dzzv[, 1], j = dzzv[, 2], x = dzzv[, 3]))
  Dzz <- matrix(nrow = dim(ridx)[1], ncol = 3)
  Dzz[, 1] <- ridx[, 1]
  Dzz[, 2] <- matrix(cidx, ncol = 1) + matrix(cs, ncol = 1)
  Dzz[, 3] <- matrix(dzzv.full, ncol = 1)
  Dzz.0 <- sparseMatrix(i = Dzz[, 1], j = Dzz[, 2], x = Dzz[, 3])
  Dzz.0 = as.matrix(Dzz.0)
  Dzz.1 <- which(Dzz.0 != 0, arr.ind = T)
  Dzz <- matrix(nrow = dim(Dzz.1)[1], ncol = 3)
  Dzz[, 1:2] <- Dzz.1[, 1:2]
  suppressMessages(Dzz[, 3] <- Dzz.0[Dzz.0!=0])
  # End Dzz
  
  # for Dxy
  dxyv.full <- as.matrix(sparseMatrix(i = dxyv[, 1], j = dxyv[, 2], x = dxyv[, 3]))
  Dxy <- matrix(nrow = dim(ridx)[1], ncol = 3)
  Dxy[, 1] <- ridx[, 1]
  Dxy[, 2] <- matrix(cidx, ncol = 1) + matrix(cs, ncol = 1)
  Dxy[, 3] <- matrix(dxyv.full, ncol = 1)
  Dxy.0 <- sparseMatrix(i = Dxy[, 1], j = Dxy[, 2], x = Dxy[, 3])
  Dxy.0 = as.matrix(Dxy.0)
  Dxy.1 <- which(Dxy.0 != 0, arr.ind = T)
  Dxy <- matrix(nrow = dim(Dxy.1)[1], ncol = 3)
  Dxy[, 1:2] <- Dxy.1[, 1:2]
  suppressMessages(Dxy[, 3] <- Dxy.0[Dxy.0!=0])
  # End Dxy
  
  # for Dxz
  dxzv.full <- as.matrix(sparseMatrix(i = dxzv[, 1], j = dxzv[, 2], x = dxzv[, 3]))
  Dxz <- matrix(nrow = dim(ridx)[1], ncol = 3)
  Dxz[, 1] <- ridx[, 1]
  Dxz[, 2] <- matrix(cidx, ncol = 1) + matrix(cs, ncol = 1)
  Dxz[, 3] <- matrix(dxzv.full, ncol = 1)
  Dxz.0 <- sparseMatrix(i = Dxz[, 1], j = Dxz[, 2], x = Dxz[, 3])
  Dxz.0 = as.matrix(Dxz.0)
  Dxz.1 <- which(Dxz.0 != 0, arr.ind = T)
  Dxz <- matrix(nrow = dim(Dxz.1)[1], ncol = 3)
  Dxz[, 1:2] <- Dxz.1[, 1:2]
  suppressMessages(Dxz[, 3] <- Dxz.0[Dxz.0!=0])
  # End Dxz
  
  # for Dyz
  dyzv.full <- as.matrix(sparseMatrix(i = dyzv[, 1], j = dyzv[, 2], x = dyzv[, 3]))
  Dyz <- matrix(nrow = dim(ridx)[1], ncol = 3)
  Dyz[, 1] <- ridx[, 1]
  Dyz[, 2] <- matrix(cidx, ncol = 1) + matrix(cs, ncol = 1)
  Dyz[, 3] <- matrix(dyzv.full, ncol = 1)
  Dyz.0 <- sparseMatrix(i = Dyz[, 1], j = Dyz[, 2], x = Dyz[, 3])
  Dyz.0 = as.matrix(Dyz.0)
  Dyz.1 <- which(Dyz.0 != 0, arr.ind = T)
  Dyz <- matrix(nrow = dim(Dyz.1)[1], ncol = 3)
  Dyz[, 1:2] <- Dyz.1[, 1:2]
  suppressMessages(Dyz[, 3] <- Dyz.0[Dyz.0!=0])
  # End Dyz

  szv <- dim(Volv)[1]
  # for Mat1gv
  diag.volv <- diag(szv)
  diag(diag.volv) <- Volv
  kron.Volv.Mat1.0 <- kron(diag.volv, Mat1)
  kron.Volv.Mat1.1 <- which(kron.Volv.Mat1.0 != 0, arr.ind = T)
  kron.Volv.Mat1.2 <- matrix(nrow = dim(kron.Volv.Mat1.1)[1], ncol = 3)
  kron.Volv.Mat1.2[, 1:2] <- kron.Volv.Mat1.1[, 1:2]
  kron.Volv.Mat1.2[, 3] <- kron.Volv.Mat1.0[kron.Volv.Mat1.0!=0]
  
  Mat1gv <- kron.Volv.Mat1.2
  # end Mat1gv

  
  # for Eg
  Mat1gv.full <- sparseMatrix(i = Mat1gv[, 1], j = Mat1gv[, 2], x = Mat1gv[, 3])
  Dxx.full.0 <- sparseMatrix(i = Dxx[, 1], j = Dxx[, 2], x = Dxx[, 3])
  Dyy.full.0 <- sparseMatrix(i = Dyy[, 1], j = Dyy[, 2], x = Dyy[, 3])
  Dzz.full.0 <- sparseMatrix(i = Dzz[, 1], j = Dzz[, 2], x = Dzz[, 3])
  Dxy.full.0 <- sparseMatrix(i = Dxy[, 1], j = Dxy[, 2], x = Dxy[, 3])
  Dxz.full.0 <- sparseMatrix(i = Dxz[, 1], j = Dxz[, 2], x = Dxz[, 3])
  Dyz.full.0 <- sparseMatrix(i = Dyz[, 1], j = Dyz[, 2], x = Dyz[, 3])
  
  Mrow <- max(dim(Dxx.full.0)[1], dim(Dyy.full.0)[1], dim(Dzz.full.0)[1], dim(Dxy.full.0)[1], dim(Dxz.full.0)[1], dim(Dyz.full.0)[1])
  Mcol <- max(dim(Dxx.full.0)[2], dim(Dyy.full.0)[2], dim(Dzz.full.0)[2], dim(Dxy.full.0)[2], dim(Dxz.full.0)[2], dim(Dyz.full.0)[2])
  
  Dxx.full <- Matrix(0, nrow = Mrow, ncol = Mcol, sparse = TRUE)
  Dyy.full <- Matrix(0, nrow = Mrow, ncol = Mcol, sparse = TRUE)
  Dzz.full <- Matrix(0, nrow = Mrow, ncol = Mcol, sparse = TRUE)
  Dxy.full <- Matrix(0, nrow = Mrow, ncol = Mcol, sparse = TRUE)
  Dxz.full <- Matrix(0, nrow = Mrow, ncol = Mcol, sparse = TRUE)
  Dyz.full <- Matrix(0, nrow = Mrow, ncol = Mcol, sparse = TRUE)
  
  Dxx.full[1:dim(Dxx.full.0)[1], 1:dim(Dxx.full.0)[2]] <- Dxx.full.0
  Dyy.full[1:dim(Dyy.full.0)[1], 1:dim(Dyy.full.0)[2]] <- Dyy.full.0
  Dzz.full[1:dim(Dzz.full.0)[1], 1:dim(Dzz.full.0)[2]] <- Dzz.full.0
  Dxy.full[1:dim(Dxy.full.0)[1], 1:dim(Dxy.full.0)[2]] <- Dxy.full.0
  Dxz.full[1:dim(Dxz.full.0)[1], 1:dim(Dxz.full.0)[2]] <- Dxz.full.0
  Dyz.full[1:dim(Dyz.full.0)[1], 1:dim(Dyz.full.0)[2]] <- Dyz.full.0
  
  Eg <- t_shallow(Dxx.full) %*% Mat1gv.full %*% Dxx.full + 
    t_shallow(Dyy.full) %*% Mat1gv.full %*% Dyy.full + 
    t_shallow(Dzz.full) %*% Mat1gv.full %*% Dzz.full + 
    2 * (t_shallow(Dxy.full) %*% Mat1gv.full %*% Dxy.full + 
           t_shallow(Dxz.full) %*% Mat1gv.full %*% Dxz.full + 
           t_shallow(Dyz.full) %*% Mat1gv.full %*% Dyz.full)
  
  return(Eg)
}
