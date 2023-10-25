thdata <- function(V, Th){
  m <- dim(Th)[1]
  if(m > 1){
    Edge <- NULL; Face <- NULL;
    TE <- NULL; TF <- NULL;
    numE <- 0; numF <- 0;
    TE.1 <- matrix(nrow = 1, ncol = 3)
    for(i in (1:m)){
      Thi <- Th[i, ]
      LocE.1 <- cbind(Thi[1], Thi[1], Thi[1], Thi[2], Thi[2], Thi[3])
      LocE.2 <- cbind(Thi[2], Thi[3], Thi[4], Thi[3], Thi[4], Thi[4])
      LocE <- rbind(LocE.1, LocE.2)
      colnames(LocE) <- c("V1", "V2", "V3", "V4", "V5", "V6")
      LocE = rbind(apply(LocE, 2, min), apply(LocE, 2, max))
      for(j in 1:6){
        if(!!length(Edge)){
          edgenum <- which((LocE[1, j] == Edge[, 1]) & (LocE[2, j] == Edge[, 2]))
        }else{
          edgenum <- {}
        }
        if(!length(edgenum)){
          Edge <- rbind(Edge, t(LocE[, j]))
          numE <- numE + 1
          edgenum <- numE
        }
        TE.1 <- matrix(nrow = 1, ncol = 3)
        TE.1[1] <- i
        TE.1[2] <- edgenum
        TE.1[3] <- 1
        TE <- rbind(TE, TE.1)
        TE <- TE[order(TE[, 2]), ]
      }
      LocF.1 <- cbind(Thi[1], Thi[2], Thi[3])
      LocF.2 <- cbind(Thi[1], Thi[2], Thi[4])
      LocF.3 <- cbind(Thi[1], Thi[3], Thi[4])
      LocF.4 <- cbind(Thi[2], Thi[3], Thi[4])
      LocF <- rbind(LocF.1, LocF.2, LocF.3, LocF.4)
      colnames(LocF) <- c("V1", "V2", "V3")
      LocF <- apply(t(LocF), 2, sort)
      for(j in 1:4){
        if (!!length(Face)){
          facenum <- which((LocF[1, j] == Face[, 1]) & (LocF[2, j] == Face[, 2]) & (LocF[3, j] == Face[, 3]))
        }else{
          facenum <- NULL
        }
        if (!length(facenum)) {
          Face <- rbind(Face, t(LocF[, j]))
          numF <- numF + 1
          facenum <- numF
        }
        TF.1 <- matrix(nrow = 1, ncol = 3)
        TF.1[1] <- i
        TF.1[2] <- facenum
        TF.1[3] <- 1
        TF <- rbind(TF, TF.1)
        TF <- TF[order(TF[, 2]), ]
      }
    }
    
    n <- dim(V)[1]
    Indx1 <- matrix(0, (4 * m), 1)
    Indx2 <- Indx1
    for (i in 1:m) {
      Indx1[(4*i-3):(4*i), 1] <- matrix(i, 4, 1)
      Indx2[(4*i-3):(4*i), 1] <- t(Th[i, ])
    }
    TV <- cbind(Indx1, Indx2, matrix(1, (4 * m), 1))
    TV <- TV[order(TV[, 2]), ]
    
    ###################################################################
    TV.full <- matrix(0, nrow = dim(Th)[1], ncol = dim(V)[1])
    TV.initial <- as.matrix(sparseMatrix(i = TV[, 1], j = TV[, 2], x = TV[, 3]))
    TV.full[1:dim(TV.initial)[1], 1:dim(TV.initial)[2]] <- TV.initial
    TV.sparse <- Matrix(TV.full, sparse = TRUE)
    ###################################################################
    
    Indx1 <- matrix(0, (2 * numE), 1)
    Indx2 <- Indx1
    for (i in 1:numE) {
      Indx1[(2*i-1):(2*i), 1] <- matrix(i, 2, 1)
      Indx2[(2*i-1):(2*i), 1] <- t(Edge[i, ])
    }
    EV <- cbind(Indx1, Indx2, matrix(1, (2 * numE), 1))
    EV.sparse <- sparseMatrix(i = EV[, 1], j = EV[, 2], x = EV[, 3])
    EV.full <- as.matrix(EV.sparse)
    EV <- EV[order(EV[, 2]), ]
    Indx1 <- matrix(0, nrow = (3 * numF), ncol = 1);
    Indx2 <- Indx1
    for (i in 1:numF) {
      Indx1[(3 * i - 2):(3 * i), 1] <- matrix(i, 3, 1)
      Indx2[(3 * i - 2):(3 * i), 1] <- t(Face[i, ])
    }
    FV <- cbind(Indx1, Indx2, matrix(1, (3 * numF), 1))
    FV.sparse <- sparseMatrix(i=FV[, 1], j = FV[, 2], x = FV[, 3])
    FV.full <- as.matrix(FV.sparse)
    FV <- FV[order(FV[, 2]), ]
    
    FE <- which(FV.full %*% t(EV.full) == 2, arr.ind = T)
    FE <- cbind(FE, matrix(1, nrow = (dim(FE)[1]), ncol = 1))
    
    BF <- matrix(nrow=(max(TF[, 2])), ncol = 3) 
    for (i in 1:(max(TF[, 2]))) {
      BF[i, 1] <- 1 
      BF[i, 2] <- i
      BF[i, 3] <- length(which(TF[, 2] == i))
    }
    BF <- which(BF[, 3] == 1) 
    
    FE.full <- sparseMatrix(i = FE[, 1], j = FE[, 2], x = FE[, 3])
    FE.full <- as.matrix(FE.full)
    dummy <- which(FE.full[BF, ] >= 1, arr.ind = TRUE)[, 1]
    BE <- which(FE.full[BF, ] >= 1, arr.ind = TRUE)[, 2]
    BE <- del_reps(BE)
    BE <- matrix(BE, ncol = 1)
    dummy <- which(EV.full[BE, ] >= 1, arr.ind = TRUE)[, 1]
    BV <- which(EV.full[BE, ] >= 1, arr.ind = TRUE)[, 2]
    BV <- del_reps(BV)
    BV <- matrix(BV, ncol = 1)
  }else{
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
                     TF = TF, TV = TV, TE = TE,
                     FV = FV, FE = FE, EV = EV,
                     BF = BF, BE = BE, BV = BV, TV.sparse = TV.sparse)
  return(thdata.list)
}