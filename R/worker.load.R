worker.load <- function(V, Tr, TV, inVT.list, Y, Z, d, nl = 1, ns, P.func = NULL) {
  
  # Automatically choose or validate the parallel backend.
  # P.func = 1: use mclapply()
  # P.func = 2: use parLapply()
  P.func <- choose.P.func(P.func)
  
  if (missing(ns) || is.null(ns) || !is.numeric(ns) || length(ns) != 1 || is.na(ns) || ns < 1) {
    stop("Argument 'ns' must be a positive integer specifying the number of cores.")
  }
  ns <- as.integer(ns)
  
  message(">>> worker.load starting with iota = ", nl)
  
  nd <- ncol(Tr)
  conditions.all <- c()
  
  if (P.func == 1) {
    
    res.all <- parallel::mclapply(
      1:nrow(Tr),
      FUN = data.sub,
      mc.cores = ns,
      V = V,
      Tr = Tr,
      TV = TV,
      inVT.list = inVT.list,
      Y = Y,
      Z = Z,
      nl = nl
    )
    
  } else if (P.func == 2) {
    
    cl <- parallel::makeCluster(ns)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = c("data.sub", "V", "Tr", "TV", "inVT.list", "Y", "Z", "nl"),
      envir = environment()
    )
    
    parallel::clusterEvalQ(cl, library(prodlim))
    
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
    nq <- (d + 2) * (d + 1) / 2
  } else if (nd == 4) {
    nq <- (d + 3) * (d + 2) * (d + 1) / 2 / 3
  } else {
    stop("'Tr' must have 3 columns for 2D triangulations or 4 columns for 3D tetrahedralizations.")
  }
  
  ind.expand <- which(conditions.all[, 1] < nq / 2)
  nl0 <- nl
  
  while (length(ind.expand) > 0) {
    
    message(">>> expanding subregions to nl0 = ", nl0 + 1)
    nl0 <- nl0 + 1
    
    if (P.func == 1) {
      
      res.all[ind.expand] <- parallel::mclapply(
        ind.expand,
        FUN = data.sub,
        mc.cores = ns,
        V = V,
        Tr = Tr,
        TV = TV,
        inVT.list = inVT.list,
        Y = Y,
        Z = Z,
        nl = nl0
      )
      
    } else if (P.func == 2) {
      
      parallel::clusterExport(
        cl,
        varlist = c("nl0"),
        envir = environment()
      )
      
      res.all[ind.expand] <- parallel::parLapply(cl, ind.expand, function(iT) {
        data.sub(iT, V, Tr, TV, inVT.list, Y, Z, nl0)
      })
    }
    
    conditions.all <- do.call("rbind", lapply(res.all, "[[", 6))
    ind.expand <- which(conditions.all[, 1] < nq / 2)
  }
  
  return(res.all)
}
