rm(list = ls())

library(devtools)
install_github("funstatpackages/MPST")
library(MPST) 
data(VT1.cube)

d = 3; r = 1; 

# Grid points for evaluation
z1.grid = seq(0.05, 0.95, by = 0.1); n1.grid = length(z1.grid); n1.grid
z2.grid = seq(0.05, 0.95, by = 0.1); n2.grid = length(z2.grid); n2.grid
z3.grid = seq(0.05, 0.95, by = 0.1); n3.grid = length(z3.grid); n3.grid
n.grid = n1.grid * n2.grid * n3.grid;
Z.grid = as.matrix(expand.grid(z1.grid, z2.grid, z3.grid))

func = 1; sigma = 1;
pop = dataGenerator3D(Z.grid, V, Tr, func, sigma)
Y.grid = pop$Y; mu.grid = pop$mu;
ind.grid = pop$ind.inside;

# Simulation parameters
n = 10000;
lambda = 10^seq(-6, 6, by = 0.5)
iter = 1
set.seed(iter)
z1 = runif(n, 0, 1)
z2 = runif(n, 0, 1)
z3 = runif(n, 0, 1)
Z = cbind(z1, z2, z3)
sam = dataGenerator3D(Z, V, Tr, func, sigma, iter)
Y = as.vector(sam$Y); Z = as.matrix(sam$Z);

###### Need to check
t0.g = proc.time()
mfit.g = fit.MPST(Y = Y, Z = Z, V = V, Tr = Tr, d = d, r = r, lambda = lambda, nl = 1, method = "G")
proc.time() - t0.g
t1.g = proc.time()
mpred.g = predict.MPST(mfit.g, Z.grid)
proc.time() - t1.g
mspe.g = mean((Y.grid - mpred.g$Ypred)^2, na.rm = TRUE)
cat("mse.g =", mfit.g$mse, "and mspe.g =", mspe.g, "\n")

t0.d = proc.time()
mfit.d = fit.MPST(Y = Y, Z = Z, V = V, Tr = Tr, d = d, r = r, lambda = lambda, nl = 1, method = "D")
proc.time() - t0.d
t1.d = proc.time()
mpred.d = predict.MPST(mfit.d, Z.grid)
proc.time() - t1.d
mspe.d = mean((Y.grid - mpred.d$Ypred)^2, na.rm = TRUE)
cat("mse.d =", mfit.d$mse, "and mspe.d =", mspe.d, "\n")

# plot.MPST(mfit.g, Z.grid)
# plot.MPST(mfit.d, Z.grid)
