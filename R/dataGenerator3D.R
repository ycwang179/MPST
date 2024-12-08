#' Data Generator for Numerical Examples in 3D
#'
#' This function generates data for numerical examples based on given 3D coordinates and a tetrahedral partition.
#'
#' @param Z A matrix of dimension \code{n} by three, where \code{n} is the number of points. Each row represents the 3D coordinates of a point.
#' @param V A matrix of dimension \code{N} by three, where \code{N} is the number of vertices in the tetrahedral mesh. Each row represents the 3D coordinates of a vertex.
#' @param Tr A matrix of dimension \code{nT} by four, where \code{nT} is the number of tetrahedrons. Each row contains the indices of four vertices in \code{V} that form a tetrahedron.
#' @param func An integer specifying the function used to calculate the mean values (\code{mu}) at the 3D points. Supported values are:
#'   \itemize{
#'     \item \code{1}: \code{(a * b + d^2) * (2 - z3^2)}
#'     \item \code{2}: \code{(a * b + d^2) * sin(pi * z3)}
#'     \item \code{3}: \code{(a * b + d^2) * cos(pi * z3)}
#'     \item \code{4}: \code{(a * b + d^2) * exp(-8 * z3^2)}
#'     \item Any other value defaults to \code{(a * b + d^2)}.
#'   }
#' @param sigma The standard deviation of the random noise added to the data. Default is \code{0.1}.
#' @param seed A seed for the random number generator to ensure reproducibility. Default is \code{2024}.
#'
#' @return A list containing the following elements:
#'   \item{Y}{A vector of generated response values, including noise.}
#'   \item{mu}{A vector of true mean values (\code{mu}) at the 3D points.}
#'   \item{Z}{The input matrix of 3D coordinates.}
#'   \item{ind.inside}{A binary vector indicating whether each point is inside (\code{1}) or outside (\code{0}) the tetrahedral partition.}
#'   \item{ind.T}{A vector containing the indices of the tetrahedron each point belongs to (if inside).}
#'
#' @examples
#' # Example with randomly generated 3D points and a simple tetrahedral partition
#' Z <- matrix(runif(300, -1, 1), ncol = 3)
#' V <- matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3, byrow = TRUE)
#' Tr <- matrix(c(1, 2, 3, 4), ncol = 4, byrow = TRUE)
#' result <- dataGenerator3D(Z, V, Tr, func = 1, sigma = 0.5, seed = 42)
#'
#' @export

dataGenerator3D <- function(Z, V, Tr, func = 1, sigma = 0.1, seed = 2024) {
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
