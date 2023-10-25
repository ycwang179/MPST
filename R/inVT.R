#' Decide whether a point is inside of a given triangulation.
#'
#' This function is used to decided whether a point is inside of a given triangulation.
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @param V The \code{nV} by two matrix of vertices of a triangulation, where \code{nV} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param Z The coordinates of 2D or 3D points. Each row is the 2D or 3D coordinates of a point.
#' \cr
#' @return A list of vectors, including:
#' \item{ind.inside}{A vector of dimension \code{n} by one matrix that lists whether the points are inside of a given triangulation. 0 -- represents outside the triangulation, while 1 -- represents inside the triangulation.}
#' \item{ind.T}{A vector contains the indexes of which triangle the points falls in.}
#' \item{lam}{A matrix contains the barycentric coordinates for each points.}
#' 
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia and Li Wang from the Iowa State University.
#'
#' @examples
#' xx=c(-0.25,0.75,0.25,1.25)
#' yy=c(-0.25,0.25,0.75,1.25)
#' V=rbind(c(0,0),c(1,0),c(1,1),c(0,1))
#' Tr=rbind(c(1,2,3),c(1,3,4))
#' inVT(V,Tr,xx,yy)
#' 
#' @export

inVT <- function(V, Tr, Z){
  Z = as.matrix(Z)
  n <- nrow(Z)
  nT <- nrow(Tr)
  tol <- -eps(1)
  nq = ncol(Z)
  if (nq == 2) {
    z1 = as.vector(Z[, 1])
    z2 = as.vector(Z[, 2])
    L = HQbary(V, Tr, z1, z2)
  } else if (nq == 3) {
    z1 = as.vector(Z[, 1])
    z2 = as.vector(Z[, 2])
    z3 = as.vector(Z[, 3])
    L = HQbary3D(V, Tr, z1, z2, z3)
  }
  
  ind.inside <- matrix(0, nrow = n, ncol = 1)
  ind.T <- matrix(NA, nrow = n, ncol = 1)
  lam <- matrix(NA, nrow = n, ncol = (nq + 1))
  for (j in 1:nT) {
    tmp = colSums(L[(j - 1) * (nq + 1) + 1:(nq + 1), ] >= tol)
    I <- which(tmp == (nq + 1))
    ind.T[I] <- j
    ind.inside[I] <- 1
    lam[I, ] <- t(L[(j - 1) * (nq + 1) + 1:(nq + 1), I])
  }
  inVT.list = list(ind.inside = ind.inside,
                   ind.T = ind.T,
                   lam = lam)
  return(inVT.list)
}