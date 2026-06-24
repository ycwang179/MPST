# MPST

**MPST** is an R package for **Multivariate Penalized Splines over Triangulations**. It provides a unified framework for smoothing, denoising, prediction, and visualization on complex two-dimensional and three-dimensional domains.

The package is designed for nonparametric regression and function estimation when observations are noisy, spatially indexed, and located on regular or irregular domains. MPST uses triangulation-based spline representations to respect complex geometry and supports both global learning for moderate-sized problems and distributed learning for larger-scale applications.

## Overview

MPST implements penalized spline smoothing over user-supplied triangulations. In two dimensions, the domain is represented by triangles; in three dimensions, it is represented by tetrahedra. The package combines:

* Bernstein polynomial spline basis construction
* Explicit smoothness constraints across simplex interfaces
* Derivative-based roughness penalization
* Global penalized spline fitting
* Distributed learning through domain decomposition
* Prediction at user-specified locations
* Built-in visualization for 2D and 3D fitted functions

This package was developed as part of my Ph.D. research in Statistics at George Mason University. The manuscript describing the package is being prepared for submission to the *Journal of Statistical Software*.

## Key Features

* **2D and 3D smoothing:** Supports bivariate and trivariate domains using triangle-based and tetrahedron-based meshes.
* **Irregular domain support:** Uses triangulations to represent non-rectangular, non-convex, and geometrically complex domains.
* **Global learning:** Fits a single penalized spline model over the full triangulated domain.
* **Distributed learning:** Partitions the domain into subregions, fits local models in parallel, and assembles the results into a globally smooth estimator.
* **Automatic smoothing parameter selection:** Supports generalized cross-validation over a user-supplied grid of smoothing parameters.
* **Prediction workflow:** Evaluates fitted models at new spatial locations or regular plotting grids.
* **Visualization tools:** Provides contour and surface plots for 2D settings and slice-based visualization for 3D settings.
* **Research software design:** Provides an end-to-end workflow for modeling, prediction, visualization, simulation studies, and scalable computation.

## When to Use MPST

MPST is useful when the data are observed over a spatial domain with complex geometry, such as:

* Irregular two-dimensional spatial regions
* Non-convex domains such as horseshoe-shaped regions
* Three-dimensional point-cloud or imaging domains
* Biomedical and medical imaging data
* Large-scale spatial or geometric datasets
* Smoothing and denoising problems where domain geometry should be respected

The package is especially helpful when regular-grid smoothers or standard tensor-product smoothers may not adequately respect the boundary or geometry of the domain.

## Installation

The development version can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ycwang179/MPST")
```

After installation, load the package with:

```r
library(MPST)
```

If the package is later released on CRAN, it can be installed with:

```r
install.packages("MPST")
```

## Basic Workflow

A typical MPST workflow consists of:

1. Preparing the response vector and spatial coordinates
2. Providing a triangulation through vertices and simplex indices
3. Fitting a global or distributed MPST model
4. Predicting fitted values at new locations
5. Visualizing the estimated function

The main user-facing functions are:

* `fit.MPST()` for model fitting
* `predict.MPST()` for prediction
* `plot.MPST()` for visualization

## Example: Global Learning

```r
library(MPST)

data(ex_square_train)

fit.g <- fit.MPST(
  y ~ m(Z, V, Tr, d = 3, r = 1),
  data = ex_square_train$m1$sigma01$tri50,
  method = "G",
  lambda = 10^seq(-6, 6, by = 0.5)
)

fit.g
```

Here, `method = "G"` fits the global MPST estimator using all observations and the full triangulation.

## Example: Prediction

```r
data(ex_square_pred)

pred.g <- predict.MPST(
  y ~ m(Z, V, Tr, d = 3, r = 1),
  data = ex_square_train$m1$sigma01$tri50,
  data.pred = ex_square_pred$m1$sigma01$tri50,
  method = "G",
  lambda = 10^seq(-6, 6, by = 0.5)
)

pred.g
```

## Example: Visualization

For two-dimensional domains, MPST supports contour and surface visualization:

```r
plot(fit.g, mview = "contour")
plot(fit.g, mview = "surface")
```

For three-dimensional domains, MPST supports slice-based visualization:

```r
plot(fit.g, mview = "slice")
```

## Example: Distributed Learning

For larger problems, MPST supports distributed learning through domain decomposition and parallel local fitting:

```r
fit.d <- fit.MPST(
  y ~ m(Z, V, Tr, d = 5, r = 1),
  data = ex_cube_train$tet48,
  method = "D",
  lambda = 10^seq(-6, 6, by = 0.5)
)

fit.d
```

In distributed learning, the domain is partitioned into subregions, local penalized spline models are fitted in parallel, and the resulting coefficients are assembled and projected to restore global smoothness.

## Main Inputs

The main inputs used by MPST are:

* `Y`: response vector
* `Z`: matrix of observation coordinates
* `V`: matrix of mesh vertices
* `Tr`: element matrix defining triangles in 2D or tetrahedra in 3D
* `d`: polynomial degree of the spline basis
* `r`: smoothness order across simplex interfaces
* `lambda`: roughness penalty parameter or candidate grid
* `method`: `"G"` for global learning or `"D"` for distributed learning

## Repository Structure

* `R/`: Main R functions for model fitting, prediction, visualization, basis construction, smoothness constraints, penalty construction, and distributed learning
* `src/`: Lower-level computational routines used by the package
* `man/`: R package documentation files
* `data/`: Example datasets used in package demonstrations
* `DESCRIPTION`: Package metadata
* `NAMESPACE`: Exported functions and package namespace configuration

## Code Samples to Review

For reviewers or employers interested in representative programming samples, I recommend reviewing:

1. **Model fitting and prediction functions** in the `R/` directory, especially functions related to `fit.MPST()`, `predict.MPST()`, and the global/distributed fitting workflows.
2. **Distributed learning components** that implement domain decomposition, local model fitting, coefficient assembly, and global smoothness projection.
3. **Basis, smoothness, and penalty construction functions** that implement the core numerical components of the MPST methodology.
4. **Computational routines in `src/`**, which demonstrate lower-level implementation for numerical efficiency.
5. **Documentation in `man/`**, which shows how the package functions are organized and documented for R users.

These components demonstrate experience with statistical programming, R package development, numerical methods, matrix computation, reproducible research software, and scalable analytical workflows.

## Technical Skills Demonstrated

This project demonstrates experience with:

* R package development
* Statistical modeling
* Nonparametric regression
* Penalized splines
* Triangulation-based methods
* Numerical optimization
* Matrix computation
* Simulation studies
* Prediction and visualization workflows
* Parallel and distributed statistical computing
* Reproducible research software development
* C/C++ integration in R package development

## Project Status

This is an active research software project developed for my Ph.D. dissertation in Statistics at George Mason University. The package and accompanying manuscript are being prepared for submission to the *Journal of Statistical Software*.

## Author

Yu-Chun Wang
Ph.D. Candidate in Statistics
George Mason University
