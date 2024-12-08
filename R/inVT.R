#' Determine whether points are inside a given triangulation.
#'
#' This function determines whether points are inside a given 2D or 3D triangulation.
#'
#' @importFrom Rcpp evalCpp
#'
#' @param V A \code{nV} by 2 or 3 matrix of vertices of the triangulation, where \code{nV} is the number of vertices. Each row contains the coordinates of a vertex.
#' @param Tr A triangulation matrix of dimension \code{nT} by 3 or 4, where \code{nT} is the number of triangles (for 2D) or tetrahedra (for 3D). Each row contains the indices of vertices in \code{V}.
#' @param Z A matrix of 2D or 3D points to check. Each row contains the coordinates of a point.
#'
#' @return A list containing:
#' \item{ind.inside}{A vector indicating whether each point is inside the triangulation. 0 represents outside, and 1 represents inside.}
#' \item{ind.T}{A vector indicating the index of the triangle (or tetrahedron) that each point falls in. NA if the point is outside.}
#' \item{lam}{A matrix containing the barycentric coordinates of each point. NA for points outside the triangulation.}
#'
#' @details This R function is adapted from the MATLAB program written by Ming-Jun Lai (University of Georgia) and Li Wang (Iowa State University).
#'
#' @examples
#' # Example: Checking if points are inside a 2D triangulation
#' xx <- c(-0.25, 0.75, 0.25, 1.25)
#' yy <- c(-0.25, 0.25, 0.75, 1.25)
#' V <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
#' Tr <- rbind(c(1, 2, 3), c(1, 3, 4))
#' Z <- cbind(xx, yy)
#' result <- inVT(V, Tr, Z)
#' print(result)
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
