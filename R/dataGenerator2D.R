#' Generate testing dataset for bivariate spline smoothing.
#'
#' @importFrom stats rnorm
#' @importFrom mgcv s
#'
#' This function generate the testing dataset for bivariate spline smoothing.
#' 
#' @param Z The cooridinates of dimension \code{n} by two. Each row is the coordinates of a point.
#' \cr
#' @param V The \code{nV} by two matrix of vertices of a triangulation, where \code{nV} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tr The triangulation matrix of dimention \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param func The choice of test function -- default is 1. Possible choices include 1, 2, 3, 4, 5, 6, 7, 8.
#' \cr
#' @param sigma The standard deviation of the white noise --  default is 0.1.
#' \cr 
#' @return A list of vectors and matrice, including:
#' \item{Y}{The response variable.}
#' \item{mu}{The mean function.}
#' \item{Z}{The coordinates.}
#' \item{ind.inside}{A vector contains the indicators whether the point is inside the given triangulation.}
#' \item{ind.T}{A vector contains the indexes of which triangle the points falls in.}
#'
#' @details This R program is modified based on Lai and Wang (2013).
#' 
#' @keywords internal


dataGenerator2D <- function(Z, V, Tr, func, sigma, seed){
  set.seed(seed)
  cat("(dataGenerator2D) sigma =", sigma, " \n")
  # location information
  if (isempty(Z)) {
    stop("The location information is required for data generation.")
  }
  Z <- matrix(Z, ncol = 2);
  z1 = as.vector(Z[, 1]); z2 = as.vector(Z[, 2]);
  n = nrow(Z)
  inVT.list <- inVT(V, Tr, Z)
  ind.inside = inVT.list$ind.inside
  ind.T = inVT.list$ind.T
  
  # test functions
  # if (!(func %in% (1:8))) {
  if (!(func %in% (1:15))) {
    stop("Test function can only be integers between 1 and 15.")
  }
  if (func == 1) {
    mu <- z1^2 + z2^2 + 2*z1*z2 # Quadratic
    hist(mu)
  }
  if (func == 2) {
    mu <- 2 * z1^3 + 2 * z2^3 # Cubic
    hist(mu)
  }
  if (func == 3) {
    mu <- 2 * z1^4 + 2 * z2^4 # Quadruplicate
    hist(mu)
  }
  if (func == 4) {
    mu <- 4 * sin(pi * z1) * sin(pi * z2) # Sine
    hist(mu)
  }
  if (func == 5) {
    mu <- (1 - cos(2*pi*z1)) * (1 - cos(2*pi*z2)) # Cosine
    hist(mu)
  }
  if (func == 6) {
    mu <- 4 * exp(-10 * ((z1-0.5)^2 + (z2-0.5)^2)) # Bump
    hist(mu)
  }
  if (func == 7) {
    mu <- 4 / (1 + exp(-10 * (z1 + z2) + 10)) # Logit
    hist(mu)
  }
  if (func == 8) {
    mu <- atan((4 * z1 - 4)^2 - (4 * z2 - 4)^2) # arctan
    hist(mu)
  }
  if (func == 9) {
    mu <- (-1) * z1^3 + z2^3 # Cubic
    hist(mu)
  }
  if (func == 10) {
    mu <- (-1) * sin(3 * pi * (z1 + 0.25)) + sin(3 * pi * z2) # Sine
    hist(mu)
    cat("function =", func, ": (-1) * sin(3 * pi * (z1 + 0.25)) + sin(3 * pi * z2) \n")
  }
  if (func == 11) {
    mu <- (-1) * sin(10 * pi * (z1 + 0.25)) + sin(10 * pi * z2) # Sine
    hist(mu)
    cat("function =", func, ": (-1) * sin(10 * pi * (z1 + 0.25)) + sin(10 * pi * z2) \n")
  }
  if (func == 12) {
    mu <- exp(-50 * ((z1-0.5)^2 + (z2-0.5)^2)) # Bump
    hist(mu)
  }
  if (func == 13) {
    mu <- 1 / (1 + exp(-10 * (z1 + z2) + 10)) # Logit
    hist(mu)
  }
  if (func == 14) {
    mu <- atan((8 * z1 - 4)^2 - (8 * z2 - 4)^2) # arctan
    hist(mu)
  }
  if (func == 15) {
    library(mgcv)
    mu <- mgcv::fs.test(z1, z2) # mgcv::fs.test(xx, yy)
    hist(mu)
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

