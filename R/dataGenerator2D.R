#' Generate Testing Dataset for Bivariate Spline Smoothing
#'
#' This function generates testing datasets for bivariate spline smoothing based on specified test functions and input coordinates.
#'
#' @importFrom mgcv fs.test
#'
#' @param Z A matrix of dimension \code{n} by two, where \code{n} is the number of points. Each row represents the 2D coordinates of a point.
#' @param V A matrix of dimension \code{nV} by two, where \code{nV} is the number of vertices in the triangulation. Each row represents the 2D coordinates of a vertex.
#' @param Tr A matrix of dimension \code{nT} by three, where \code{nT} is the number of triangles in the triangulation. Each row contains the indices of three vertices in \code{V} that form a triangle.
#' @param func An integer specifying the function used to calculate the mean values (\code{mu}) at the 2D points. Supported values range from \code{1} to \code{15}, with different mathematical forms for each choice. Default is \code{1}.
#' @param sigma The standard deviation of the Gaussian noise added to the response variable. Default is \code{0.1}.
#' @param seed A seed for the random number generator to ensure reproducibility. Default is \code{2024}.
#'
#' @return A list containing the following elements:
#'   \item{Y}{A vector of generated response values, including noise.}
#'   \item{mu}{A vector of true mean values (\code{mu}) at the 2D points.}
#'   \item{Z}{The input matrix of 2D coordinates.}
#'   \item{ind.inside}{A binary vector indicating whether each point is inside (\code{1}) or outside (\code{0}) the triangulation.}
#'   \item{ind.T}{A vector containing the indices of the triangle each point belongs to (if inside).}
#'
#' @details This R program is based on methods described in Lai and Wang (2013). It supports up to 15 test functions, each providing different mathematical forms for the mean function (\code{mu}).
#'
#' @examples
#' # Example with randomly generated 2D points and a simple triangulation
#' Z <- matrix(runif(200, -1, 1), ncol = 2)
#' V <- matrix(c(0, 0, 1, 0, 0, 1), ncol = 2, byrow = TRUE)
#' Tr <- matrix(c(1, 2, 3), ncol = 3, byrow = TRUE)
#' result <- dataGenerator2D(Z, V, Tr, func = 1, sigma = 0.1, seed = 42)
#'
#' @export

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
    # hist(mu)
  }
  if (func == 2) {
    mu <- 2 * z1^3 + 2 * z2^3 # Cubic
    # hist(mu)
  }
  if (func == 3) {
    mu <- 2 * z1^4 + 2 * z2^4 # Quadruplicate
    # hist(mu)
  }
  if (func == 4) {
    mu <- 4 * sin(pi * z1) * sin(pi * z2) # Sine
    # hist(mu)
  }
  if (func == 5) {
    mu <- (1 - cos(2*pi*z1)) * (1 - cos(2*pi*z2)) # Cosine
    # hist(mu)
  }
  if (func == 6) {
    mu <- 4 * exp(-10 * ((z1-0.5)^2 + (z2-0.5)^2)) # Bump
    # hist(mu)
  }
  if (func == 7) {
    mu <- 4 / (1 + exp(-10 * (z1 + z2) + 10)) # Logit
    # hist(mu)
  }
  if (func == 8) {
    mu <- atan((4 * z1 - 4)^2 - (4 * z2 - 4)^2) # arctan
    # hist(mu)
  }
  if (func == 9) {
    mu <- (-1) * z1^3 + z2^3 # Cubic
    # hist(mu)
  }
  if (func == 10) {
    mu <- (-1) * sin(3 * pi * (z1 + 0.25)) + sin(3 * pi * z2) # Sine
    # hist(mu)
    # cat("function =", func, ": (-1) * sin(3 * pi * (z1 + 0.25)) + sin(3 * pi * z2) \n")
  }
  if (func == 11) {
    mu <- (-1) * sin(10 * pi * (z1 + 0.25)) + sin(10 * pi * z2) # Sine
    # hist(mu)
    # cat("function =", func, ": (-1) * sin(10 * pi * (z1 + 0.25)) + sin(10 * pi * z2) \n")
  }
  if (func == 12) {
    mu <- exp(-50 * ((z1-0.5)^2 + (z2-0.5)^2)) # Bump
    # hist(mu)
  }
  if (func == 13) {
    mu <- 1 / (1 + exp(-10 * (z1 + z2) + 10)) # Logit
    # hist(mu)
  }
  if (func == 14) {
    mu <- atan((8 * z1 - 4)^2 - (8 * z2 - 4)^2) # arctan
    # hist(mu)
  }
  if (func == 15) {
    mu <- mgcv::fs.test(z1, z2)
    # hist(mu)
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

