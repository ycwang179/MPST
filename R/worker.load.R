worker.load <- function (V, Tr, TV, inVT.list, Y, Z, d, nl = 1, ns = 1) {
  nd = ncol(Tr)
  
  conditions.all <- c()
  load.all <- vector("list", nrow(Tr))

  res.all <- parallel::mclapply(1:nrow(Tr), FUN = data.sub, mc.cores = ns, 
                                V = V, Tr = Tr, TV = TV, 
                                inVT.list = inVT.list, Y = Y, Z = Z, nl = nl)
  
  # for (iT in 1:nrow(Tr)) {
  #   load.all[[iT]] = data.sub(iT, V, Tr, TV, inVT.list, Y, Z, nl = 1)
  # }
  
  conditions.all <- do.call("rbind", lapply(res.all, "[[", 6))
  
  if (nd == 3) {
    nq = (d + 2) * (d + 1) / 2
  } else if (nd == 4) {
    nq = (d + 3) * (d + 2) * (d + 1) / 2 / 3
  }
  
  ind.expand <- which(conditions.all[, 1] < nq)
  nl0 <- nl
  while (length(ind.expand) > 0) {
    nl0 = nl0 + 1
    res.all[ind.expand] <- parallel::mclapply(ind.expand, FUN = data.sub, mc.cores = ns, 
                                              V, Tr, TV, inVT.list, Y, Z, nl0)
    conditions.all <- do.call("rbind", lapply(res.all, "[[", 8))
    ind.expand <- which(conditions.all[, 1] < nq)
  }

  return(res.all)
}