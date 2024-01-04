#' Bivariate Spline Basis Function
#'
#' This function generates the basis for bivariate spline over triangulation.
#'
#' @importFrom Matrix Matrix
#' @importFrom pracma isempty
#' @importFrom Rcpp evalCpp
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 5, and usually \code{d} is greater than one. -1 represents piecewise constant.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param Z The cooridinates of dimension \code{n} by two. Each row is the coordinates of a point.
#' \cr
#' @param Hmtx The indicator of whether the smoothness matrix \code{H} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param Kmtx The indicator of whether the energy matrix \code{K} need to be generated -- default is \code{TRUE}.
#' \cr
#' @param QR The indicator of whether a QR decomposition need to be performed on the smoothness matrix -- default is \code{TRUE}.
#' \cr
#' @param TA The indicator of whether the area of the triangles need to be calculated -- default is \code{TRUE}.
#' \cr
#' @return A list of vectors and matrice, including:
#' \item{B}{The spline basis function of dimension \code{n} by \code{nT}*\code{{(d+1)(d+2)/2}}, where \code{n} is the number of observationed points, \code{nT} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.}
#' \item{ind.inside}{A vector contains the indexes of all the points which are inside the triangulation.}
#' \item{H}{The smoothness matrix.}
#' \item{Q2}{The Q2 matrix after QR decomposition of the smoothness matrix \code{H}.}
#' \item{K}{The thin-plate energy function.}
#' \item{tria.all}{The area of each triangle within the given triangulation.}
#'
#' @details This R program is modified based on the Matlab program written by Ming-Jun Lai from the University of Georgia and Li Wang from the Iowa State University.
#'
#' @examples
#' # example 1
#' xx = c(-0.25, 0.75, 0.25, 1.25)
#' yy = c(-0.25, 0.25, 0.75,1 .25)
#' Z = cbind(xx, yy)
#' d = 4; r = 1;
#' V0 = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
#' Tr0 = rbind(c(1, 2, 3), c(1, 3, 4))
#' basis(V0, Tr0, d, r, Z)
#' 
#' @export

basis.d <- function(V, Tr, d = 5, r = 1, Z){
  V <- as.matrix(V); dim(V)
  Tr <- as.matrix(Tr); dim(Tr)
  nT = nrow(Tr)
  nq = choose(d + 2, 2)
  Z <- matrix(Z, ncol = 2)
  nZ <- nrow(Z)
  
  inVT.list = inVT(V, Tr, Z)
  ind.T = inVT.list$ind.T
  Lam = inVT.list$lam
  
  # sz = sum(inVT.list$ind.inside == 1);
  # tmp = HQgetInd(V, Tr, Z[, 1], Z[, 2])
  # tmp1 = rep(1, nT)
  # tmp2 = diag(rep(1, nq))
  # c = kronecker(tmp1, tmp2)
  # IndP = 1:sz
  # A = matrix(0, sz, nq)
  # for (ii in 1:nq) {
  #   tmp4 = c[, ii]
  #   tmp3=seval(V, Tr, tmp4, sz, IndP, tmp$cnt, tmp$Lam);
  #   A[, ii] = tmp3;
  # }
  # B.test = HQblkdiag(A, tmp$cnt); dim(B.test)
  # Bi <- BSpline(V, Tr, d, r, Z[, 1], Z[, 2])$Bi; dim(Bi)
  # sum(B.test == Bi)

  B = Matrix(0, nrow = nZ, ncol = nT * nq, sparse = TRUE)
  for (i in sort(unique(ind.T))) {
    Tr.i = matrix(Tr[i, ], ncol = 3); dim(Tr.i)
    V.i = matrix(V[Tr.i, ], ncol = 2); dim(V.i)
    ind.i = (1:nZ)[ind.T == i]
    Z.i = matrix(Z[ind.i, ], ncol = 2); dim(Z.i)
    nZ.i = length(ind.i)
    Lam.i = Lam[ind.i, ]; dim(Lam.i)
    
    Ind1 <- rep(0, choose(d + 2, 3))
    Ind2 <- Ind1
    Ind3 <- Ind1
    cnt1 <- rep(0, d + 1)
    for (l in 1:d) {
      cnt1[l + 1] <- cnt1[l] + (d - l + 1) * (d - l + 2) / 2
      IJK <- indices(d - l + 1)
      I <- IJK$I
      J <- IJK$J
      K <- IJK$K
      IJK1 <- indices(d - l)
      I1 <- IJK1$I
      J1 <- IJK1$J
      K1 <- IJK1$K
      Ind1[(cnt1[l] + 1):cnt1[l + 1]] <- locate(I1 + 1, J1, K1, I, J, K)
      Ind2[(cnt1[l] + 1):cnt1[l + 1]] <- locate(I1, J1 + 1, K1, I, J, K)
      Ind3[(cnt1[l] + 1):cnt1[l + 1]] <- locate(I1, J1, K1 + 1, I, J, K)
    }
    
    A = matrix(0, nZ.i, nq)
    if (nZ.i > 0) {
      for (j in 1:nq) {
        tmp1 = rep(0, nq)
        tmp1[j] = 1
        tmp2 = seval.d(tmp1, Ind1, Ind2, Ind3, cnt1, Lam.i) 
        A[, j] = tmp2;
      }
    }
    
    B[ind.i, ((i - 1) * nq + 1:nq)] = A
  }
  
  ind.inside = which(inVT.list$ind.inside == 1) 
  basis.list = list(B = B, ind.inside = ind.inside, ind.T = ind.T)
  return(basis.list)
}
