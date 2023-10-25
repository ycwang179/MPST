data.sub <- function(iT, V, Tr, TV, inVT.list, Y, Z, nl = 1) {
  cat("Working on triangle:", iT, "\n")
  
  # obtain subregions
  Tr0 = Tr[iT, ]
  VT.list = findNabr(Tr0, V, Tr, TV, nl)
  Vs = VT.list$V1; Trs = VT.list$Tr1; t0 = VT.list$Tr0
  # Triangulation::TriPlot(Vs, Trs)
  # Triangulation:::tplot(V, Tr0, col = "blue", lwd = 3)
  # Triangulation:::tplot(V, VT.list$Tr.sub[1, ], col = 2, lwd = 2)
  # Triangulation:::tplot(V, VT.list$Tr.sub[2, ], col = 3, lwd = 2)
  # Triangulation:::tplot(V, VT.list$Tr.sub[3, ], col = 4, lwd = 2)
  # Triangulation:::tplot(V, VT.list$Tr.sub[4, ], col = 5, lwd = 2)
  # Triangulation:::tplot(V, VT.list$Tr.sub[5, ], col = 6, lwd = 2)
  # Triangulation:::tplot(V, VT.list$Tr.sub[6, ], col = 7, lwd = 2)
  # Triangulation:::tplot(V, VT.list$Tr.sub[7, ], col = 8, lwd = 2)
  # Triangulation:::tplot(V, t0, lwd = 2)

  # obtain data within each sub-domain
  Tr.sub = VT.list$Tr.sub
  ind.T1 = prodlim::row.match(as.data.frame(Tr.sub), as.data.frame(Tr))
  # inVT.list = inVT(V, Tr, Z)
  ind.T2 = inVT.list$ind.T
  n = nrow(Z)
  inds = (1:n)[ind.T2 %in% ind.T1]
  Ys = Y[inds]
  Zs = Z[inds, ] 
  # inds0 = (1:n)[iT %in% ind.T1]
  # check whether expand the sub-region.
  # C1: whether there are enough observation points for the sub-domain
  # C2: check number of triangles within each sub-domain
  conditions <- c(nrow(Zs)/nrow(Trs), nrow(Trs))

  data.sub.list = list(Vs = Vs, 
                       Trs = Trs, 
                       t0 = t0, 
                       ind.T = ind.T1,
                       inds = inds,
                       # Ys = Ys,
                       # Zs = Zs, 
                       conditions = conditions)
  return(data.sub.list)
}
