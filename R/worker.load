worker.load <- function (V, Tr, TV, inVT.list, Y, Z, d, nl = 1, ns, P.func) {
  nd = ncol(Tr)
  
  conditions.all <- c()
  
  if (P.func == 1) {
    res.all <- parallel::mclapply(1:nrow(Tr), FUN = data.sub, mc.cores = ns, 
                                  V = V, Tr = Tr, TV = TV, 
                                  inVT.list = inVT.list, Y = Y, Z = Z, nl = nl)
  } else if (P.func == 2) {
    cl <- parallel::makeCluster(ns)
    parallel::clusterExport(cl, varlist = c("data.sub", "V", "Tr", "TV", "inVT.list", "Y", "Z", "nl"), envir = environment())
    parallel::clusterEvalQ(cl, library(prodlim))
    
    conditions.all <- c()
    res.all <- parallel::parLapply(cl, 1:nrow(Tr), function(iT) {
      data.sub(iT, V, Tr, TV, inVT.list, Y, Z, nl)
    })
  }

  # load.all <- vector("list", nrow(Tr))
  # for (iT in 1:nrow(Tr)) {
  #   load.all[[iT]] = data.sub(iT, V, Tr, TV, inVT.list, Y, Z, nl = 1)
  # }
  # res.all = load.all
  conditions.all <- do.call("rbind", lapply(res.all, "[[", 6))
  
  if (nd == 3) {
    nq = (d + 2) * (d + 1) / 2
  } else if (nd == 4) {
    nq = (d + 3) * (d + 2) * (d + 1) / 2 / 3
  }
  
  ind.expand <- which(conditions.all[, 1] < nq/2)
  nl0 <- nl
  while (length(ind.expand) > 0) {
    nl0 = nl0 + 1
    if (P.func == 1){
      res.all[ind.expand] <- parallel::mclapply(ind.expand, FUN = data.sub, mc.cores = ns, 
                                                V, Tr, TV, inVT.list, Y, Z, nl0)
    } else if (P.func == 2) {
      parallel::clusterExport(cl, varlist = c("nl0"), envir = environment())
      res.all[ind.expand] <- parallel::parLapply(cl, ind.expand, function(iT) {
        data.sub(iT, V, Tr, TV, inVT.list, Y, Z, nl0)
      })
    }
    
    conditions.all <- do.call("rbind", lapply(res.all, "[[", 6))
    ind.expand <- which(conditions.all[, 1] < nq/2)
  }
  
  if (P.func == 2) {
    parallel::stopCluster(cl)
  }
  return(res.all)
}
