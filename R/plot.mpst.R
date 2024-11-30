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
  } else if ((mview == "contour") && (nd == 4)) {
    fig <- plot.contour.mpst(x, Zgrid = Zgrid, ...)
    print(fig)
  } else if (mview == "surface") {
    fig <- plot.surface.mpst(x, Zgrid = Zgrid, ...)
    print(fig)
  } else if ((mview == "slice") && (nd == 4)) {
    fig <- plot.slice.mpst(x, Zgrid = Zgrid, ...)
    print(fig)
  } else {
    stop("Invalid mview or unsupported dimensionality")
  }
  invisible(x)
}

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
      return(list(Zgrid = Zgrid, u1 = z1.grid, v1 = z2.grid, v2 = z3.grid))
    } else {
      stop("Unsupported dimensionality for grid initialization.")
    }
  } else {
    # Use the existing Zgrid to extract coordinates
    if (nd == 3) {
      u1 <- sort(unique(Zgrid[, 1]))
      v1 <- sort(unique(Zgrid[, 2]))
      return(list(Zgrid = Zgrid, u1 = u1, v1 = v1))
    } else if (nd == 4) {
      u1 <- sort(unique(Zgrid[, 1]))
      v1 <- sort(unique(Zgrid[, 2]))
      v2 <- sort(unique(Zgrid[, 3]))
      return(list(Zgrid = Zgrid, u1 = u1, v1 = v1, v2 = v2))
    } else {
      stop("Unsupported dimensionality for grid initialization.")
    }
  }
}

plot.contour.mpst <- function(mfit, Zgrid = NULL) {
  if (!("Tr" %in% names(mfit)) || !("Z" %in% names(mfit))) {
    stop("The input object 'mfit' must contain 'Tr' and 'Z'.")
  }
  
  grid_info <- initialize.grid(mfit, Zgrid)
  Zgrid <- grid_info$Zgrid
  u1 <- grid_info$u1
  v1 <- grid_info$v1
  
  mpred <- pred.MPST(mfit, Znew = Zgrid)
  if (!("Ypred" %in% names(mpred))) {
    stop("'pred.MPST()' did not return the expected 'Ypred' component.")
  }
  
  z1 <- matrix(mpred$Ypred, length(u1), length(v1), byrow = TRUE)
  
  fig <- plot_ly(
    type = "contour",
    x = u1, 
    y = v1, 
    z = t(z1), 
    contours = list(showlabels = TRUE)
  )
  return(fig)
}

plot.surface.mpst <- function(mfit, Zgrid = NULL) {
  if (!("Tr" %in% names(mfit)) || !("Z" %in% names(mfit))) {
    stop("The input object 'mfit' must contain 'Tr' and 'Z'.")
  }
  
  grid_info <- initialize.grid(mfit, Zgrid)
  Zgrid <- grid_info$Zgrid
  u1 <- grid_info$u1
  v1 <- grid_info$v1
  
  mpred <- pred.MPST(mfit, Znew = Zgrid)
  if (!("Ypred" %in% names(mpred))) {
    stop("'pred.MPST()' did not return the expected 'Ypred' component.")
  }
  
  z1 <- matrix(mpred$Ypred, nrow = length(u1), byrow = TRUE)
  df1 <- data.frame(x = mfit$Z[, 1], y = mfit$Z[, 2], z = mfit$Y)
  
  fig <- plot_ly(
    x = u1, 
    y = v1, 
    z = t(z1),
    type = "surface"
  ) %>% 
    add_trace(
      data = df1,
      x = ~x, y = ~y, z = ~z,
      mode = "markers", type = "scatter3d",
      marker = list(size = 1, color = "black")
    ) %>%
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

plot.slice.mpst <- function(mfit, Zgrid = NULL) {
  if (!("Tr" %in% names(mfit)) || !("Z" %in% names(mfit))) {
    stop("The input object 'mfit' must contain 'Tr' and 'Z'.")
  }
  
  grid_info <- initialize.grid(mfit, Zgrid, n1 = 181, n2 = 91, n3 = 46)
  Zgrid <- grid_info$Zgrid
  z1.grid <- grid_info$u1
  z2.grid <- grid_info$v1
  z3.grid <- grid_info$v2
  
  mpred <- pred.MPST(mfit, Znew = Zgrid)
  if (!("Ypred" %in% names(mpred))) {
    stop("'pred.MPST()' did not return the expected 'Ypred' component.")
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
    print(paste("Color choice:", color_choice))  
    base_palette <- switch(color_choice,
                           "1" = gray.colors(64, start = 0.2, end = 0.9),
                           "2" = rainbow(64),
                           "3" = heat.colors(64),
                           "4" = terrain.colors(64),
                           "5" = topo.colors(64),
                           "6" = cm.colors(64))
    print(base_palette[1:3])  
    return(base_palette)
  }
  
  # Define the function to plot slices
  plot_slices <- function(axial_slice, coronal_slice, sagittal_slice, color = 4) {
    col_palette <- select_color_palette(color)  
    par(mfrow = c(1, 3))  
    
    # Draw the axial slice
    image.plot(1:dim.size[1], 1:dim.size[2], new.array[,,axial_slice], 
               main = paste("Axial Plane (z =", axial_slice, ")"), 
               xlab = "x", ylab = "y", col = col_palette, axes = FALSE, useRaster = TRUE)
    
    # Draw the coronal slice
    image.plot(1:dim.size[1], 1:dim.size[3], new.array[,coronal_slice,], 
               main = paste("Coronal Plane (y =", coronal_slice, ")"), 
               xlab = "x", ylab = "z", col = col_palette, axes = FALSE, useRaster = TRUE)
    
    # Draw the sagittal slice
    image.plot(1:dim.size[2], 1:dim.size[3], new.array[sagittal_slice,,], 
               main = paste("Sagittal Plane (x =", sagittal_slice, ")"), 
               xlab = "y", ylab = "z", col = col_palette, axes = FALSE, useRaster = TRUE)
  }
  
  # Use manipulate to create an interactive interface
  manipulate(
    plot_slices(axial_slice, coronal_slice, sagittal_slice, color),
    axial_slice = slider(1, dim.size[3], initial = dim.size[3] %/% 2, label = "Axial Slice"),
    coronal_slice = slider(1, dim.size[2], initial = dim.size[2] %/% 2, label = "Coronal Slice"),
    sagittal_slice = slider(1, dim.size[1], initial = dim.size[1] %/% 2, label = "Sagittal Slice"),
    color = slider(1, 6, initial = 1, label = "Color Table (1-Gray, 2-Rainbow, 3-Heat, 4-Terrain, 5-Topo, 6-Cyan Magenta)")
  )
}
