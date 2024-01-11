rm(list = ls())

library(MPST) 
data(VT1.hs3D)

d = 4; r = 1; 

# Grid points for evaluation
z1.grid = seq(-0.85, 3.4, by = 0.1); n1.grid = length(z1.grid); n1.grid
z2.grid = seq(-0.85, 0.85, by = 0.1); n2.grid = length(z2.grid); n2.grid
z3.grid = seq(-0.5, 0.5, by = 0.1); n3.grid = length(z3.grid); n3.grid
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

t0 = proc.time()
B.original <- basis3D(V, Tr, d, r, Z)$B
t1 = proc.time() - t0
t1
cat("Time for generating the spline basis via original method is", t1[3], "\n")

t0 = proc.time()
B <- basis3D.d(V, Tr, d, r, Z)$B
t1 = proc.time() - t0
t1
cat("Time for generating the spline basis via distributed method is", t1[3], "\n")

dim(B.original); dim(B);
max(abs(B.original - B))

t0 = proc.time()
H <- smoothness3D(V, Tr, d, r)
t1 = proc.time() - t0
t1
cat("Time for calculate the smootheness constraints is", t1[3], "\n")

t0 = proc.time()
K <- energy3D(V, Tr, d)
t1 = proc.time() - t0
t1
cat("Time for generating the energy function is", t1[3], "\n")

