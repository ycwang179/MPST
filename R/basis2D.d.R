#' Basis Function Generator for Bivariate Spline over Triangulation
#'
#' @importFrom Matrix sparseMatrix det
#' @importFrom pracma eps
#'
#' @param V The \code{nV} by two matrix of vertices of a triangle, where \code{nV} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangular partition matrix of dimension \code{nT} by three, where \code{nT} is the number of triangles in the partition. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 9, and usually \code{d} is greater than one. -1 represents piecewise constant.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param Z The coordinates of dimension \code{n} by two. Each row is the 2D coordinates of a point.
#' \cr
#' @return A list of vectors and matrices, including:
#' \item{B}{The spline basis function of dimension \code{n} by \code{nT}*\code{{(d+2)(d+1)/2}}, where \code{n} is the number of observed points, \code{nT} is the number of triangles in the given partition, and \code{d} is the degree of polynomials for the spline. If some points do not fall in the triangular partition, the generation of the spline basis will not take those points into consideration.}
#' \item{ind.inside}{A vector contains the indexes of all the points which are inside the triangulation.}
#' \item{ind.T}{A vector contains the indexes of the triangles each point falls in.}
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia.
#'
#' @examples
#' # example 1
#' d = 2; r = 1;
#' x1.o = seq(0.1, 0.9, by = 0.1)
#' x2.o = seq(0.1, 0.9, by = 0.1)
#' Z = expand.grid(x1.o, x2.o)
#' V <- matrix(rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1)), ncol = 2)
#' Tr <- matrix(rbind(c(1, 2, 3), c(2, 3, 4)), ncol = 3)
#' B.all = basis2D.d(V, Tr, d, r, Z)
#'
#' # example 2
#' d = 3; r = 1;
#' x1.o = seq(0.1, 0.9, by = 0.2)
#' x2.o = seq(0.1, 0.9, by = 0.2)
#' Z = expand.grid(x1.o, x2.o)
#' V <- matrix(rbind(c(0, 0), c(0, 1), c(1, 0), c(1, 1), c(0.5, 0.5)), ncol = 2)
#' Tr <- matrix(rbind(c(1, 2, 5), c(2, 4, 5), c(4, 3, 5), c(3, 1, 5)), ncol = 3)
#' B.all = basis2D.d(V, Tr, d, r, Z)
#' @export
basis2D.d <- function (V, Tr, d, r, Z) {
  
  Z <- as.matrix(Z);
  nq <- choose(d + 2, 2)
  nT <- nrow(Tr)
  
  inVT.list = inVT(V, Tr, Z)
  ind.inside = which(inVT.list$ind.inside == 1)
  ind.T = inVT.list$ind.T
  lam = inVT.list$lam
  
  ind.T = ind.T[ind.inside]
  lam = lam[ind.inside, ]
  sind.T <- sort(ind.T)
  Z = Z[ind.inside, ]
  nZ = nrow(Z)
  
  B = Matrix(0, nrow = nZ, ncol = nT * nq, sparse = TRUE)
  tmp = expand.grid(d:0, 0:d, 0:d)
  exps <- as.matrix(tmp[rowSums(tmp) == d, ]) 
  ne = nrow(exps)
  div <- factorial(exps[, 1]) * factorial(exps[, 2]) * factorial(exps[, 3])
  
  for (i in sort(unique(ind.T))) {
    Tr.i = matrix(Tr[i, ], ncol = 3); dim(Tr.i)
    V.i = matrix(V[Tr.i, ], ncol = 2); dim(V.i)
    ind.i = (1:nZ)[ind.T == i]
    Z.i = matrix(Z[ind.i, ], ncol = 2); dim(Z.i)
    nZ.i = length(ind.i)
  
    val1.i <- sweep(t(matrix(rep(lam[ind.i, 1], ne), ncol = ne)), 1, exps[, 1], "^")
    val2.i <- sweep(t(matrix(rep(lam[ind.i, 2], ne), ncol = ne)), 1, exps[, 2], "^")
    val3.i <- sweep(t(matrix(rep(lam[ind.i, 3], ne), ncol = ne)), 1, exps[, 3], "^")
    val4.i <- val1.i * val2.i * val3.i
    val5.i <- sweep(val4.i, 1, div, "/")
    val.i <- val5.i * factorial(d)
    
    B[ind.i, ((i - 1) * nq + 1:nq)] = t(val.i)
  }
  
  basis.list <- list(B = B,
                     ind.inside = ind.inside,
                     ind.T = ind.T,
                     inVT.list = inVT.list)
  
  return(basis.list)
}
