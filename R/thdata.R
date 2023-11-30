thdata <- function(V, Tr){
  nT <- nrow(Tr)
  if(nT > 1){
    Edge <- NULL; Face <- NULL;
    TE.ind <- NULL; TF.ind <- NULL;
    numE <- 0; numF <- 0;
    
    for (i in 1:nT) {
      Tri <- Tr[i, ]
      LocE <- rbind(c(Tri[1], Tri[1], Tri[1], Tri[2], Tri[2], Tri[3]), 
                    c(Tri[2], Tri[3], Tri[4], Tri[3], Tri[4], Tri[4]))
      LocE = rbind(apply(LocE, 2, min), apply(LocE, 2, max))
      
      for (j in 1:6) {
        if (!is.null(Edge)) {
          edgenum <- which((LocE[1, j] == Edge[, 1]) & (LocE[2, j] == Edge[, 2]))
        } else {
          edgenum <- NULL
        }
        if (length(edgenum) == 0) {
          Edge <- rbind(Edge, t(LocE[, j]))
          numE <- numE + 1
          edgenum <- numE
        }
        TE.ind <- rbind(TE.ind, c(i, edgenum, 1))
        # TE.ind <- TE.ind[order(TE.ind[, 2]), ]
      }
      
      LocF <- rbind(c(Tri[1], Tri[2], Tri[3]), 
                    c(Tri[1], Tri[2], Tri[4]), 
                    c(Tri[1], Tri[3], Tri[4]), 
                    c(Tri[2], Tri[3], Tri[4]))
      LocF <- apply(t(LocF), 2, sort)
      for (j in 1:4) {
        if (!is.null(Face)) {
          facenum <- which((LocF[1, j] == Face[, 1]) & (LocF[2, j] == Face[, 2]) & (LocF[3, j] == Face[, 3]))
        } else {
          facenum <- NULL
        }
        if (length(facenum) == 0) {
          Face <- rbind(Face, t(LocF[, j]))
          numF <- numF + 1
          facenum <- numF
        }
        TF.ind <- rbind(TF.ind, c(i, facenum, 1))
        # TF.ind <- TF.ind[order(TF.ind[, 2]), ]
      }
    }
    TE = sparseMatrix(i = TE.ind[, 1], j = TE.ind[, 2], x = TE.ind[, 3])
    TF = sparseMatrix(i = TF.ind[, 1], j = TF.ind[, 2], x = TF.ind[, 3])
    
    nV <- nrow(V)
    Indx1 <- matrix(0, (4 * nT), 1)
    Indx2 <- Indx1
    for (i in 1:nT) {
      Indx1[(4*i-3):(4*i), 1] <- matrix(i, 4, 1)
      Indx2[(4*i-3):(4*i), 1] <- t(Tr[i, ])
    }
    TV.ind <- cbind(Indx1, Indx2, matrix(1, (4 * nT), 1))
    TV = sparseMatrix(i = TV.ind[, 1], j = TV.ind[, 2], x = TV.ind[, 3])
    
    Indx1 <- matrix(0, (2 * numE), 1)
    Indx2 <- Indx1
    for (i in 1:numE) {
      Indx1[(2*i-1):(2*i), 1] <- matrix(i, 2, 1)
      Indx2[(2*i-1):(2*i), 1] <- t(Edge[i, ])
    }
    EV.ind <- cbind(Indx1, Indx2, matrix(1, (2 * numE), 1))
    EV <- sparseMatrix(i = EV.ind[, 1], j = EV.ind[, 2], x = EV.ind[, 3])
    
    Indx1 <- matrix(0, nrow = (3 * numF), ncol = 1);
    Indx2 <- Indx1
    for (i in 1:numF) {
      Indx1[(3 * i - 2):(3 * i), 1] <- matrix(i, 3, 1)
      Indx2[(3 * i - 2):(3 * i), 1] <- t(Face[i, ])
    }
    FV.ind <- cbind(Indx1, Indx2, matrix(1, (3 * numF), 1))
    FV <- sparseMatrix(i = FV.ind[, 1], j = FV.ind[, 2], x = FV.ind[, 3])
    
    tmp = FV %*% t_shallow(EV)
    ind.tmp = which(tmp@x == 2)
    j.tmp = rep(1:tmp@Dim[2], diff(tmp@p))
    FE.ind <- cbind(tmp@i[ind.tmp] + 1, j.tmp[ind.tmp])
    FE.ind <- cbind(FE.ind, matrix(1, nrow = nrow(FE.ind), ncol = 1))
    FE <- sparseMatrix(i = FE.ind[, 1], j = FE.ind[, 2], x = FE.ind[, 3])
    
    BF = which(Matrix::colSums(TF) == 1)
    # BF <- matrix(nrow = (max(TF.ind[, 2])), ncol = 3) 
    # for (i in 1:(max(TF.ind[, 2]))) {
    #   BF[i, 1] <- 1 
    #   BF[i, 2] <- i
    #   BF[i, 3] <- length(which(TF.ind[, 2] == i))
    # }
    # BF <- which(BF[, 3] == 1) 
    
    BE <- which(as.matrix(FE)[BF, ] > 0, arr.ind = TRUE)[, 2]
    BE = unique(BE)
    # BE <- del_reps(BE)
    
    BV <- which(as.matrix(EV)[BE, ] > 0, arr.ind = TRUE)[, 2]
    BV = unique(BV)
    # BV <- del_reps(BV)
  } else {
    Edge <- matrix(c(1, 1, 1, 2, 2, 3, 2, 3, 4, 3, 4, 4), nrow = 6, ncol = 2)
    Face <- matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4), nrow = 4, ncol = 3)
    TF <- matrix(1, nrow = 1, ncol = 4)
    TE <- matrix(1, nrow = 1, ncol = 6)
    TV <- matrix(1, nrow = 1, ncol = 4)
    FV <- matrix(1, nrow = 4, ncol = 4)
    FE <- matrix(1, nrow = 4, ncol = 6)
    EV <- matrix(c(1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1), nrow = 6, ncol = 4)
    BF <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
    BE <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 6, ncol = 1)
    BV <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  }
  
  thdata.list = list(Edge = Edge, Face = Face,
                     TF.ind = TF.ind, TF = TF,
                     TV.ind = TV.ind, TV = TV,
                     TE.ind = TE.ind, TE = TE,
                     FV.ind = FV.ind, FV = FV, 
                     FE.ind = FE.ind, FE = FE,
                     EV.ind = EV.ind, EV = EV,
                     BF = BF, BE = BE, BV = BV)
  return(thdata.list)
}
