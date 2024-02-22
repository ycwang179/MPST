rm(list = ls())

# library(devtools)
# install_github("funstatpackages/MPST")
library(MPST)  
library(profvis)
data(VT.square)
# Triangulation::TriPlot(V, Tr)

d = 5; r = 1; 

# Grid points for evaluation
n1.grid = 51; n2.grid = 51; n.grid = n1.grid * n2.grid;
u.grid = seq(0, 1, length.out = n1.grid)
v.grid = seq(0, 1, length.out = n2.grid)
uu.grid = rep(u.grid, each = n2.grid)
vv.grid = rep(v.grid, times = n1.grid)
Z.grid = as.matrix(cbind(uu.grid, vv.grid))

func = 1; sigma = 1;
pop = dataGenerator2D(Z.grid, V, Tr, func, sigma)
Y.grid = pop$Y; mu.grid = pop$mu;
ind.grid = pop$ind.inside;

# Simulation parameters
n = 5000;
lambda = 10^seq(-6, 6, by = 0.5)
iter = 1
Z = matrix(runif(2*n, 0, 1), nrow = n, ncol = 2)
sam = dataGenerator2D(Z, V, Tr, func, sigma, iter)
Y = as.vector(sam$Y); Z = as.matrix(sam$Z);

t0.g = proc.time()
profvis({
mfit.g = fit.MPST(Y = Y, Z = Z, V = V, Tr = Tr, d = d, r = r, lambda = lambda, method = "G")
})
proc.time() - t0.g

t1.g = proc.time()
mpred.g = predict.MPST(mfit.g, Z.grid)
proc.time() - t1.g
mspe.g = mean((Y.grid - mpred.g$Ypred)^2, na.rm = TRUE)
cat("mse.g =", mfit.g$mse, "and mspe.g =", mspe.g, "\n")

t0.d = proc.time()
profvis({
mfit.d = fit.MPST(Y = Y, Z = Z, V = V, Tr = Tr, d = d, r = r, lambda = lambda, nl = 1, method = "D")
})
proc.time() - t0.d
mpred.d = predict.MPST(mfit.d, Z.grid)
proc.time() - t1.d
mspe.d = mean((Y.grid - mpred.d$Ypred)^2, na.rm = TRUE)
cat("mse.d =", mfit.d$mse, "and mspe.d =", mspe.d, "\n")

# plot.MPST(mfit.g, Z.grid)
# plot.MPST(mfit.d, Z.grid)
