#' Basis Function for Trivariate Spline over Triangulation
#'
#' @importFrom Matrix sparseMatrix det
#' @importFrom pracma eps
#'
#' @param V The \code{nV} by three matrix of vertices of a tetrahedron, where \code{nV} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The tetrahedral partition matrix of dimension \code{nT} by four, where \code{nT} is the number of tetrahedrons in the partition. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 9, and usually \code{d} is greater than one. -1 represents piecewise constant.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param Z The coordinates of dimension \code{n} by three. Each row is the 3D coordinates of a point.
#' \cr
#' @return A list of vectors and matrices, including:
#' \item{B}{The spline basis function of dimension \code{n} by \code{nT}*\code{{(d+3)(d+2)(d+1)/(2*3)}}, where \code{n} is the number of observed points, \code{nT} is the number of tetrahedrons in the given partition, and \code{d} is the degree of polynomials for the spline. If some points do not fall in the tetrahedral partition, the generation of the spline basis will not take those points into consideration.}
#' \item{ind.inside}{A vector contains the indexes of all the points which are inside the triangulation.}
#' \item{ind.T}{A vector contains the indexes of the triangles each point falls in.}
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia.
#'
#' @examples
#' # example 1
#' d = 2; r = 1;
#' z1.o = seq(0.1, 0.9, by = 0.1)
#' z2.o = seq(0.1, 0.9, by = 0.1)
#' z3.o = seq(0.1, 0.9, by = 0.1)
#' Z = expand.grid(z1.o, z2.o, z3.o)
#' V <- matrix(rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 0), c(1, 1, 0),
#' c(0, 0, 1), c(0, 1, 1), c(1, 0, 1), c(1, 1, 1)), ncol = 3)
#' Th <- matrix(rbind(c(5, 1, 2, 3), c(6, 5, 2, 3), c(6, 7, 5, 3),
#' c(6, 4, 7, 3), c(6, 2, 4, 3), c(6, 8, 7, 4)), ncol = 4)
#' B.all = basis3D(V, Th, d, r, Z)
#'
#' # example 2
#' d = 2; r = 1;
#' z1.o = seq(0.1, 0.9, by = 0.2)
#' z2.o = seq(0.1, 0.9, by = 0.2)
#' z3.o = seq(0.1, 0.9, by = 0.2)
#' Z = expand.grid(z1.o, z2.o, z3.o)
#' V <- matrix(rbind(c(0, 0, 0), c(0, 0.5, 0), c(0, 1, 0),
#' c(0.5, 0, 0), c(0.5, 0.5, 0), c(0.5, 1, 0),
#' c(1, 0, 0), c(1, 0.5, 0), c(1, 1, 0),
#' c(0, 0, 0.5), c(0, 0.5, 0.5), c(0, 1, 0.5),
#' c(0.5, 0, 0.5), c(0.5, 0.5, 0.5), c(0.5, 1, 0.5),
#' c(1, 0, 0.5), c(1, 0.5, 0.5), c(1, 1, 0.5),
#' c(0, 0, 1), c(0, 0.5, 1), c(0, 1, 1),
#' c(0.5, 0, 1), c(0.5, 0.5, 1), c(0.5, 1, 1),
#' c(1, 0, 1), c(1, 0.5, 1), c(1, 1, 1)), ncol = 3)
#' Th <- matrix(rbind(c(16, 22, 23, 25), c(15, 17, 23, 24),
#' c(3, 5, 6, 12), c(14, 15, 21, 23), c(12, 14, 15, 21),
#' c(11, 12, 14, 20), c(5, 11, 12, 14), c(14, 20, 21, 23),
#' c(14, 20, 22, 23), c(6, 12, 14, 15), c(5, 6, 12, 14),
#' c(18, 24, 26, 27), c(12, 14, 20, 21), c(7, 13, 14, 16),
#' c(15, 21, 23, 24), c(14, 16, 17, 23), c(2, 4, 5, 11),
#' c(2, 3, 5, 11), c(3, 5, 11, 12), c(5, 7, 13, 14),
#' c(4, 5, 7, 13), c(17, 18, 24, 26), c(17, 23, 25, 26),
#' c(16, 17, 23, 25), c(14, 15, 17, 23), c(14, 16, 22, 23),
#' c(13, 14, 20, 22), c(13, 14, 16, 22), c(17, 23, 24, 26),
#' c(2, 4, 10, 11), c(1, 2, 4, 10), c(6, 8, 14, 15),
#' c(4, 5, 11, 13), c(5, 7, 8, 14), c(13, 19, 20, 22),
#' c(10, 11, 13, 19), c(4, 10, 11, 13), c(5, 6, 8, 14),
#' c(9, 15, 17, 18), c(15, 17, 18, 24), c(5, 11, 13, 14),
#' c(11, 13, 14, 20), c(11, 13, 19, 20), c(8, 14, 16, 17),
#' c(7, 8, 14, 16), c(8, 14, 15, 17), c(8, 9, 15, 17),
#' c(6, 8, 9, 15)), ncol = 4)
#' B.all = basis3D(V, Th, d, r, Z)
#' @export

basis3D <- function (V, Tr, d, r, Z) {
  # memory.size(max = TRUE) # Windows-specific function
  Z <- as.matrix(Z);
  m <- choose(d + 3, 3)
  Tr <- t(apply(t(Tr), 2, sort))
  nT <- nrow(Tr)

  inVT.list = inVT(V, Tr, Z)
  ind.inside = inVT.list$ind.inside
  ind.T = inVT.list$ind.T
  lam = inVT.list$lam

  # pointlocation3D.list <- pointLocation3D(Z, V, Tr)
  # ind.T <- pointlocation3D.list$ind.VT
  # lam <- pointlocation3D.list$lam
  # ind.inside <- which(ind.T != 0)

  Z <- Z[ind.inside == 1, ]; dim(Z)
  n <- nrow(Z)
  lam <- lam[ind.inside == 1, ]
  ind.T <- matrix(ind.T[ind.inside == 1], ncol = 1)
  sind.T <- sort(ind.T)
  ind <- matrix(order(ind.T), ncol = 1)

  b1 <- matrix(lam[ind, 1], ncol = 1)
  b2 <- matrix(lam[ind, 2], ncol = 1)
  b3 <- matrix(lam[ind, 3], ncol = 1)
  b4 <- matrix(lam[ind, 4], ncol = 1)

  exps <- loop(d); ne = nrow(exps);
  div <- matrix(factorial(exps[, 1]) * factorial(exps[, 2]) * factorial(exps[, 3]) * factorial(exps[, 4]), ncol = 1)
  exps.1 <- matrix(exps[, 1], ncol = 1)
  exps.2 <- matrix(exps[, 2], ncol = 1)
  exps.3 <- matrix(exps[, 3], ncol = 1)
  exps.4 <- matrix(exps[, 4], ncol = 1)

  val.1 <- sweep(t(matrix(rep(t(b1), ne), ncol = ne)), 1, exps.1, "^")
  val.2 <- sweep(t(matrix(rep(t(b2), ne), ncol = ne)), 1, exps.2, "^")
  val.3 <- sweep(t(matrix(rep(t(b3), ne), ncol = ne)), 1, exps.3, "^")
  val.4 <- sweep(t(matrix(rep(t(b4), ne), ncol = ne)), 1, exps.4, "^")
  val.5 <- val.1 * val.2 * val.3 * val.4
  val.6 <- sweep(val.5, 1, div, "/")
  val <- val.6 * factorial(d)

  rind1 <- matrix(rep(matrix((1:m), nrow = 1), n), ncol = n)
  row.shift <- m * (t(sind.T) - 1)
  rind <- sweep(rind1, 2, row.shift, "+")
  cind <- t(matrix(rep((1:n), m), nrow = n))
  Bi.full <- matrix(0, m * nT, n)
  Bi <- as.matrix(sparseMatrix(i = as.vector(rind), j = as.vector(cind), x = as.vector(val)))
  Bi.full[1:dim(Bi)[1], 1:dim(Bi)[2]] <- Bi
  B <- matrix(0, m * nT, n)
  B[, ind] <- as.matrix(Bi.full)
  Bi <- t(as.matrix(Bi)); B <- t(B);
  Bi <- Matrix(Bi, sparse = TRUE)
  B <- Matrix(B, sparse = TRUE)

  basis.list <- list(B = B,
                     Bi = Bi,
                     ind = ind,
                     ind.inside = ind.inside,
                     ind.T = ind.T)

  return(basis.list)
}
