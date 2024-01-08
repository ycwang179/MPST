rm(list = ls())

# library(devtools)
# install_github("funstatpackages/MPST")
library(MPST)  

# install_github("funstatpackages/Triangulation")
library(Triangulation)

# library(profvis)
data(VT.square)
Triangulation::TriPlot(V, Tr)

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
  
t0 = proc.time()
inVT.list = inVT(V, Tr, Z)
t1 = proc.time() - t0
cat("Time for checking the point locations relative to the triangulation is", t1[3], "\n")
  
ind.inside = which(inVT.list$ind.inside == 1)
ind.nna <- (1:n)[!is.na(Y)]
ind = intersect(ind.inside, ind.nna)
Zi <- Z[ind, ]
Yi <- Y[ind]
ni = length(Yi)
  
t0 = proc.time()
B.original <- basis(V, Tr, d, r, Zi)$B
t1 = proc.time() - t0
cat("Time for generating the spline basis via original method is", t1[3], "\n")

t0 = proc.time()
B <- basis.d(V, Tr, d, r, Zi)$B
t1 = proc.time() - t0
cat("Time for generating the spline basis via distributed method is", t1[3], "\n")

t0 = proc.time()
B.new <- basis2D.d(V, Tr, d, r, Zi)$B
t1 = proc.time() - t0
cat("Time for generating the spline basis via new distributed method is", t1[3], "\n")

dim(B.original); dim(B); dim(B.new)
max(abs(B.original - B))
max(abs(B.original - B.new))
