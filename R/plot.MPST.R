#' Plot Different Views for MPST Models
#'
#' @importFrom plotly plot_ly
#' @importFrom manipulate manipulate
#' @importFrom fields image.plot
#' @importFrom magrittr %>%
#'
#' @description This function allows plotting different views (contour, surface, or slice) 
#' for an MPST model fit object.
#' @param x An MPST model fit object.
#' @param Zgrid An optional grid for plotting. If NULL, it will be generated automatically.
#' @param mview A character string specifying the view to plot. Must be one of 'contour', 
#' 'surface', or 'slice'.
#' @param ... Additional arguments passed to the specific plot functions.
#' @return No return value. The function generates plots.
#' @examples
#' \dontrun{
#' plot.MPST(model_fit_object, mview = "contour")
#' }
#' @export
plot.MPST <- function(x, Zgrid = NULL, mview = NULL, ...) {
  if (!inherits(x, "MPST")) {
    stop("Object must be of class 'MPST'")
  }
  if (is.null(mview)) {
    stop("mview must be specified (e.g., 'contour', 'surface', or 'slice').")
  }
  if (!mview %in% c("contour", "surface", "slice")) {
    stop("Invalid mview. Must be one of 'contour', 'surface', or 'slice'.")
  }
   
  nd <- ncol(x$Tr)
  if ((mview == "contour") && (nd == 3)) {
    fig <- plot.contour.mpst(x, Zgrid = Zgrid, ...)
    print(fig)
  } else if ((mview == "surface") && (nd == 3)) {
    fig <- plot.surface.mpst(x, Zgrid = Zgrid, ...)
    print(fig)
  } else if ((mview == "slice") && (nd == 4)) {
    plot.slice.mpst(x, Zgrid = Zgrid, ...)
  } else {
    stop("Invalid mview or unsupported dimensionality")
  }
  invisible()
}

#' Initialize Grid for MPST Plotting
#'
#' @description Generates a grid for plotting based on the input MPST model and dimensionality.
#' This function is used internally and is not intended for direct use by end users.
#' @param mfit An MPST model fit object.
#' @param Zgrid An optional grid. If NULL, a grid is created automatically.
#' @param n1 Number of points in the first dimension.
#' @param n2 Number of points in the second dimension.
#' @param n3 Number of points in the third dimension (if applicable).
#' @return A list containing the generated grid and corresponding coordinate vectors.
#' @keywords internal
initialize.grid <- function(mfit, Zgrid = NULL, n1 = 101, n2 = 101, n3 = 101) {
  nd <- ncol(mfit$Tr)
  if (is.null(Zgrid)) {
    # Calculate the range of each dimension
    limits <- apply(mfit$Z[, 1:(nd-1), drop = FALSE], 2, function(col) {
      c(min(col) - 0.0001, max(col) + 0.0001)
    })
    
    if (nd == 3) {
      z1.grid <- seq(limits[1, 1], limits[2, 1], length.out = n1)
      z2.grid <- seq(limits[1, 2], limits[2, 2], length.out = n2)
      Zgrid <- as.matrix(expand.grid(z1.grid, z2.grid))
      return(list(Zgrid = Zgrid, u1 = z1.grid, v1 = z2.grid))
    } else if (nd == 4) {
      z1.grid <- seq(limits[1, 1], limits[2, 1], length.out = n1)
      z2.grid <- seq(limits[1, 2], limits[2, 2], length.out = n2)
      z3.grid <- seq(limits[1, 3], limits[2, 3], length.out = n3)
      Zgrid <- as.matrix(expand.grid(z1.grid, z2.grid, z3.grid))
      return(list(Zgrid = Zgrid, u1 = z1.grid, v1 = z2.grid, w1 = z3.grid))
    } else {
      stop("Unsupported dimensionality for grid initialization.")
    }
  } else {
    # Use the existing Zgrid to extract coordinates
    if (nd == 3) {
      z1.grid <- sort(unique(Zgrid[, 1]))
      z2.grid <- sort(unique(Zgrid[, 2]))
      return(list(Zgrid = Zgrid, u1 = z1.grid, v1 = z2.grid))
    } else if (nd == 4) {
      z1.grid <- sort(unique(Zgrid[, 1]))
      z2.grid <- sort(unique(Zgrid[, 2]))
      z3.grid <- sort(unique(Zgrid[, 3]))
      return(list(Zgrid = Zgrid, u1 = z1.grid, v1 = z2.grid, w1 = z3.grid))
    } else {
      stop("Unsupported dimensionality for grid initialization.")
    }
  }
}

#' Plot Contour for MPST Models
#'
#' @description Generates a contour plot for a 3D MPST model. This function is used internally
#' by `plot.MPST` and is not intended for direct use by end users.
#' @param mfit An MPST model fit object.
#' @param Zgrid An optional grid for plotting. If NULL, a grid will be generated.
#' @return A plotly object representing the contour plot.
#' @keywords internal
plot.contour.mpst <- function(mfit, Zgrid = NULL) {
  if (!("Tr" %in% names(mfit)) || !("Z" %in% names(mfit))) {
    stop("The input object 'mfit' must contain 'Tr' and 'Z'.")
  }
  
  grid_info <- initialize.grid(mfit, Zgrid)
  Zgrid <- grid_info$Zgrid
  u1 <- grid_info$u1
  v1 <- grid_info$v1
  
  mpred <- pred.mpst(mfit, Znew = Zgrid)
  if (!("Ypred" %in% names(mpred))) {
    stop("'pred.mpst()' did not return the expected 'Ypred' component.")
  }
  
  z1 <- matrix(mpred$Ypred, length(u1), length(v1), byrow = TRUE)
  
  fig <- plotly::plot_ly(
    type = "contour",
    x = u1, 
    y = v1, 
    z = t(z1), 
    contours = list(showlabels = TRUE)
  )
  return(fig)
}

#' Plot Surface for MPST Models
#'
#' @description Generates a 3D surface plot for an MPST model. This function is used internally
#' by `plot.MPST` and is not intended for direct use by end users.
#' @param mfit An MPST model fit object.
#' @param Zgrid An optional grid for plotting. If NULL, a grid is generated automatically.
#' @return A plotly object representing the surface plot.
#' @keywords internal
plot.surface.mpst <- function(mfit, Zgrid = NULL) {
  if (!("Tr" %in% names(mfit)) || !("Z" %in% names(mfit))) {
    stop("The input object 'mfit' must contain 'Tr' and 'Z'.")
  }
  
  grid_info <- initialize.grid(mfit, Zgrid)
  Zgrid <- grid_info$Zgrid
  u1 <- grid_info$u1
  v1 <- grid_info$v1
  
  mpred <- pred.mpst(mfit, Znew = Zgrid)
  if (!("Ypred" %in% names(mpred))) {
    stop("'pred.mpst()' did not return the expected 'Ypred' component.")
  }
  
  z1 <- matrix(mpred$Ypred, nrow = length(u1), byrow = TRUE)
  df1 <- data.frame(x = mfit$Z[, 1], y = mfit$Z[, 2], z = mfit$Y)
  
  fig <- plotly::plot_ly(
    x = u1, 
    y = v1, 
    z = t(z1),
    type = "surface"
  ) magrittr::%>% 
    add_trace(
      data = df1,
      x = ~x, y = ~y, z = ~z,
      mode = "markers", type = "scatter3d",
      marker = list(size = 1, color = "black")
    ) magrittr::%>%
    layout(
      scene = list(
        camera = list(eye = list(x = -1.6, y = -1.6, z = 0.8)),
        xaxis = list(title = "Z1"),
        yaxis = list(title = "Z2"),
        zaxis = list(title = "Value")
      )
    )
  
  return(fig)
}

#' Plot Slices for 3D MPST Models
#'
#' @description Generates interactive 2D slices for a 3D MPST model. This function is used internally
#' by `plot.MPST` and is not intended for direct use by end users.
#' @param mfit An MPST model fit object.
#' @param Zgrid An optional grid for plotting. If NULL, a grid is generated automatically.
#' @return A manipulate object that provides interactive slices through the 3D array.
#' @keywords internal
plot.slice.mpst <- function(mfit, Zgrid = NULL) {
  if (!("Tr" %in% names(mfit)) || !("Z" %in% names(mfit))) {
    stop("The input object 'mfit' must contain 'Tr' and 'Z'.")
  }
  
  grid_info <- initialize.grid(mfit, Zgrid, n1 = 181, n2 = 91, n3 = 46)
  Zgrid <- grid_info$Zgrid
  z1.grid <- grid_info$u1
  z2.grid <- grid_info$v1
  z3.grid <- grid_info$w1
  
  if (anyNA(Zgrid)) {
    stop("Generated Zgrid contains NA values. Please check input data.")
  }
  
  mpred <- pred.mpst(mfit, Znew = Zgrid)
  if (!("Ypred" %in% names(mpred))) {
    stop("'pred.mpst()' did not return the expected 'Ypred' component.")
  }
  
  # Fill 3D array
  indices.matrix <- expand.grid(i = 1:length(z1.grid), j = 1:length(z2.grid), k = 1:length(z3.grid))
  indices.matrix <- as.matrix(indices.matrix)
  new.array <- array(NA, dim = c(length(z1.grid), length(z2.grid), length(z3.grid)))
  
  for (row in 1:nrow(indices.matrix)) {
    coords <- indices.matrix[row, ]
    new.array[coords[1], coords[2], coords[3]] <- mpred$Ypred[row]
  }
  
  # Interactive slices
  dim.size <- dim(new.array)
  
  # Define the function to select a color palette
  select_color_palette <- function(color_choice) {
    base_palette <- switch(color_choice,
                           "1" = gray.colors(64),
                           "2" = rainbow(64),
                           "3" = heat.colors(64),
                           "4" = terrain.colors(64),
                           "5" = topo.colors(64),
                           "6" = cm.colors(64))
    return(base_palette)
  }
  
  # Define the function to plot slices
  plot_slices <- function(axial_slice, coronal_slice, sagittal_slice, color = 2) {
    col_palette <- select_color_palette(color)
    par(mfrow = c(1, 3))
    
    # Axial slice
    if (all(is.na(new.array[, , axial_slice]))) {
      warning("Axial slice contains only NA values.")
      plot.new()
      text(0.5, 0.5, "Axial slice contains only NA values.", cex = 1.5)
    } else {
      fields::image.plot(1 : dim.size[1], 1 : dim.size[2], new.array[, , axial_slice], 
                 main = paste("Axial Plane (z =", axial_slice, ")"), 
                 xlab = "x", ylab = "y", col = col_palette, axes = FALSE, useRaster = TRUE)
    }
    
    # Coronal slice
    if (all(is.na(new.array[, coronal_slice, ]))) {
      warning("Coronal slice contains only NA values.")
      plot.new()
      text(0.5, 0.5, "Coronal slice contains only NA values.", cex = 1.5)
    } else {
      fields::image.plot(1 : dim.size[1], 1 : dim.size[3], new.array[, coronal_slice, ], 
                 main = paste("Coronal Plane (y =", coronal_slice, ")"), 
                 xlab = "x", ylab = "z", col = col_palette, axes = FALSE, useRaster = TRUE)
    }
    
    # Sagittal slice
    if (all(is.na(new.array[sagittal_slice, , ]))) {
      warning("Sagittal slice contains only NA values.")
      plot.new()
      text(0.5, 0.5, "Sagittal slice contains only NA values.", cex = 1.5)
    } else {
      fields::image.plot(1 : dim.size[2], 1 : dim.size[3], new.array[sagittal_slice, , ], 
                 main = paste("Sagittal Plane (x =", sagittal_slice, ")"), 
                 xlab = "y", ylab = "z", col = col_palette, axes = FALSE, useRaster = TRUE)
    }
  }
  
  # Use manipulate to create an interactive interface
  manipulate::manipulate(
    plot_slices(axial_slice, coronal_slice, sagittal_slice, color),
    axial_slice = slider(1, max(1, dim.size[3]), initial = max(1, dim.size[3] %/% 2), label = "Axial Slice"),
    coronal_slice = slider(1, max(1, dim.size[2]), initial = max(1, dim.size[2] %/% 2), label = "Coronal Slice"),
    sagittal_slice = slider(1, max(1, dim.size[1]), initial = max(1, dim.size[1] %/% 2), label = "Sagittal Slice"),
    color = slider(1, 6, initial = 1, label = "Color Table (1-Gray, 2-Rainbow, 3-Heat, 4-Terrain, 5-Topo, 6-Cyan Magenta)")
  )
  invisible(NULL)
}
