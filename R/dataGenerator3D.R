#' Data Generator for Numerical Example
#'
#' @param V The \code{N} by three matrix of vertices of a tetrahedron, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The tetrahedral partition matrix of dimension \code{nT} by four, where \code{nT} is the number of tetrahedrons in the partition. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param Z The coordinates of dimension \code{n} by three. Each row is the 3D coordinates of a point.
#' \cr
#' @return A data frame.
#'
#' @export
dataGenerator3D <- function(Z, V, Tr, func = 1, sigma = 1, seed = 2023) {
  set.seed(seed)
  
  # location information
  if (isempty(Z)) {
    stop("The location information is required for data generation.")
  }
  
  Z <- matrix(Z, ncol = 3);
  z1 = as.vector(Z[, 1]); z2 = as.vector(Z[, 2]); z3 = as.vector(Z[, 3]);
  n = nrow(Z)
  
  inVT.list <- inVT(V, Tr, Z)
  ind.inside = inVT.list$ind.inside
  ind.T = inVT.list$ind.T
  
  r0 = 0.1; r = 0.5; l = 3; b = 1;
  q = pi * r/2;

  mu = rep(NA, n);
  a = rep(0, n); d = rep(0, n);

  # Part 1
  ind = which((z1 >= 0) & (z2 > 0));
  a[ind] = q + z1[ind];
  d[ind] = z2[ind] - r;

  # Part 2
  ind = which((z1 >= 0) & (z2 <= 0));
  a[ind] = -q - z1[ind];
  d[ind] = -r - z2[ind];

  # Part 3
  ind = which(z1 < 0);
  a[ind] = -atan(z2[ind]/z1[ind]) * r;
  d[ind] = sqrt(z1[ind]^2 + z2[ind]^2) - r;

  # ind = which((abs(d) > r - r0) | (z1 > l & (z1 - l)^2 + d^2 > (r - r0)^2));
  # ind = (1:n)[ind.inside == 1]
  
  if (func == 1) {
    mu = (a * b + d^2) * (2 - z3^2);
  } else if (func == 2) {
    mu = (a * b + d^2) * sin(pi * z3);
  } else if (func == 3) {
    mu = (a * b + d^2) * cos(pi * z3);
  } else if (func == 4) {
    mu = (a * b + d^2) * exp(-8 * z3^2);
  } else {
    mu = (a * b + d^2);
  }

  eps <- rnorm(n, mean = 0, sd = sigma)
  mu[ind.inside == 0] <- NA; eps[ind.inside == 0] <- NA;
  Y <- mu + eps
  
  dat = list(Y = Y, 
             mu = mu, 
             Z = Z, 
             ind.inside = ind.inside, 
             ind.T  = ind.T)
  return(dat)
}
