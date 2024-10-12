rm(list = ls())

#library(devtools)
#install_github("funstatpackages/MPST")
library(MPST)  
data(VT1.hs3D)
# Triangulation::TriPlot(V, Tr)

d = 3; r = 1; 

# Grid points for evaluation
z1.grid = seq(-0.85, 3.4, by = 0.05); n1.grid = length(z1.grid); n1.grid
z2.grid = seq(-0.85, 0.85, by = 0.05); n2.grid = length(z2.grid); n2.grid
z3.grid = seq(-0.5, 0.5, by = 0.05); n3.grid = length(z3.grid); n3.grid
n.grid = n1.grid * n2.grid * n3.grid;
Z.grid = as.matrix(expand.grid(z1.grid, z2.grid, z3.grid))

func = 1; sigma = 1;
pop = dataGenerator3D(Z.grid, V, Tr, func, sigma)
Y.grid = pop$Y; mu.grid = pop$mu;
ind.grid = pop$ind.inside;

# Simulation parameters
n = 50000;
lambda = 10^seq(-6, 6, by = 0.5)
iter = 1
z1 = runif(2 * n, -0.8, 3.35)
z2 = runif(2 * n, -0.8, 0.8)
z3 = runif(2 * n, -0.45, 0.45)
Z = cbind(z1, z2, z3)

t0 = proc.time()
inVT.list = inVT(V, Tr, Z)
t1 = proc.time() - t0
cat("Time for checking the point locations relative to the triangulation is", t1[3], "\n")

Z = Z[inVT.list$ind.inside == 1, ]
Z = Z[1:n, ]
sam = dataGenerator3D(Z, V, Tr, func, sigma, iter)
Y = as.vector(sam$Y); Z = as.matrix(sam$Z);

t0.g = proc.time()
mfit.g = fit.MPST(Y = Y, Z = Z, V = V, Tr = Tr, d = d, r = r, lambda = lambda, nl = 1, method = "G")
proc.time() - t0.g
mpred.g = predict.MPST(mfit.g, Z.grid)
mspe.g = mean((Y.grid - mpred.g$Ypred)^2, na.rm = TRUE)
cat("mse.g =", mfit.g$mse, "and mspe.g =", mspe.g, "\n")

t0.d = proc.time()
mfit.d = fit.MPST(Y = Y, Z = Z, V = V, Tr = Tr, d = d, r = r, lambda = lambda, nl = 1, method = "D")
proc.time() - t0.d
mpred.d = predict.MPST(mfit.d, Z.grid)
mspe.d = mean((Y.grid - mpred.d$Ypred)^2, na.rm = TRUE)
cat("mse.d =", mfit.d$mse, "and mspe.d =", mspe.d, "\n")
