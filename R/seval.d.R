seval.d <- function(cval, Ind1, Ind2, Ind3, cnt1, Lam.i){
  Lam1 = Lam.i[, 1]
  Lam2 = Lam.i[, 2]
  Lam3 = Lam.i[, 3]
  C <- cval[Ind1[(cnt1[1] + 1):cnt1[2]]] %*% t(Lam1) +
    cval[Ind2[(cnt1[1] + 1):cnt1[2]]] %*% t(Lam2) +
    cval[Ind3[(cnt1[1] + 1):cnt1[2]]] %*% t(Lam3)
  
  s1 <- Matrix(diag(Lam1), sparse = TRUE); dim(s1)
  s2 <- Matrix(diag(Lam2), sparse = TRUE); dim(s2)
  s3 <- Matrix(diag(Lam3), sparse = TRUE); dim(s3)
  
  if (d >= 2) {
    for (k in 2:d) {
      C <- C[Ind1[(cnt1[k] + 1):cnt1[k + 1]], ] %*% s1 +
        C[Ind2[(cnt1[k] + 1):cnt1[k + 1]], ] %*% s2 +
        C[Ind3[(cnt1[k] + 1):cnt1[k + 1]], ] %*% s3
    }
  }
  return(as.vector(C))
}
