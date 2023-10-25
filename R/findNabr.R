findNabr <- function(Tr0, V, Tr, TV, nl = 1){
  nd = ncol(Tr) 
  # TV = tdata(V, Tr)
  # TV <- as.matrix(TV$TV)
  # V.sub = matrix(Tr0, ncol = 1)
  V.sub = c(Tr0)
  V1 = c(); Tr1 = c();
  
  if (nl == 0) {
    V1 = V[Tr0, ]; 
    Tr1 = matrix(1:nd, ncol = nd);
    V.sub = prodlim::row.match(Tr0, Tr);
    Tr.sub = matrix(Tr0, ncol = nd);
    Tr0 = Tr1; 
  } else { 
    for (i in 1:nl){
      k = length(V.sub)
      J = c()
      for(j in 1:k){
        J1 = which(TV[, V.sub[j]] != 0)
        J = c(J, J1)
        # J1 = matrix(which(TV[, V.sub[j, 1]] != 0, arr.ind = TRUE), ncol = 1)
        # J <- rbind(J, J1)
      }
      V2 = Tr[J, ]
      V.sub = c(V2)
      V1 = c(V1, V.sub)
      Tr1 = c(Tr1, J)
    }
    V1 = sort(unique(V1))
    V.sub = V1
    V1 = V[V.sub, ]
    Tr1 = Tr[sort(unique(Tr1)), ]
    Tr.sub = Tr1
    nT1 = nrow(Tr1)
    for(i in 1:nT1){
      if (nd == 3) {
        j = Tr1[i, 1]; nj = which(V.sub == j); Tr1[i, 1] = nj
        j = Tr1[i, 2]; nj = which(V.sub == j); Tr1[i, 2] = nj
        j = Tr1[i, 3]; nj = which(V.sub == j); Tr1[i, 3] = nj
      } else if (nd == 4) {
        j = Tr1[i, 1]; nj = which(V.sub == j); Tr1[i, 1] = nj
        j = Tr1[i, 2]; nj = which(V.sub == j); Tr1[i, 2] = nj
        j = Tr1[i, 3]; nj = which(V.sub == j); Tr1[i, 3] = nj
        j = Tr1[i, 4]; nj = which(V.sub == j); Tr1[i, 4] = nj
      }
    }
    if (nd == 3) {
      j = Tr0[1]; nj = which(V.sub == j); Tr0[1] = nj
      j = Tr0[2]; nj = which(V.sub == j); Tr0[2] = nj
      j = Tr0[3]; nj = which(V.sub == j); Tr0[3] = nj
    } else if (nd == 4) {
      j = Tr0[1]; nj = which(V.sub == j); Tr0[1] = nj
      j = Tr0[2]; nj = which(V.sub == j); Tr0[2] = nj
      j = Tr0[3]; nj = which(V.sub == j); Tr0[3] = nj
      j = Tr0[4]; nj = which(V.sub == j); Tr0[4] = nj
    }
  }
  
  findNabr.list = list(V1 = V1, 
                       Tr1 = Tr1, 
                       Tr0 = Tr0, 
                       V.sub = V.sub, 
                       Tr.sub = Tr.sub)
  return(findNabr.list)
}