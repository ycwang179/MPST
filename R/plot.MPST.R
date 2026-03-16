#' Plot Contour, Surface, or Slice Views for MPST Models
#'
#' @description
#' Generates contour, surface, or slice visualizations from a fitted
#' \pkg{MPST} model object. Rather than directly plotting fitted values at the
#' original training locations, this function evaluates the fitted model on a
#' plotting grid and visualizes the resulting values. If \code{Zgrid} is not
#' supplied, a plotting grid is generated automatically from the domain of the
#' fitted data.
#'
#' For two-dimensional domains, \code{mview = "contour"} and
#' \code{mview = "surface"} produce contour and surface plots, respectively.
#' For three-dimensional domains, \code{mview = "slice"} produces orthogonal
#' slice views through the fitted volume.
#'
#' In addition, precomputed 3D arrays can be plotted directly for slice
#' visualization when the input object contains \code{Yarray},
#' \code{xgrid}, \code{ygrid}, and \code{zgrid}.
#'
#' @param x
#' An object of class \code{"MPST"}, typically returned by
#' \code{\link{fit.MPST}}. For slice plotting, it may also be an MPST-like
#' object containing precomputed components \code{Yarray}, \code{xgrid},
#' \code{ygrid}, and \code{zgrid}.
#' @param Zgrid
#' An optional matrix of plotting locations. If \code{NULL}, a plotting grid is
#' generated automatically from the coordinate range of the fitted data.
#' @param mview
#' A character string specifying the type of visualization. Must be one of
#' \code{"contour"}, \code{"surface"}, or \code{"slice"}.
#' @param slice_style
#' A character string specifying the slice plotting style. Supported values are
#' \code{"default"} and \code{"brain"}. This argument is only used when
#' \code{mview = "slice"}.
#' @param ...
#' Additional arguments passed to the internal plotting routines.
#'
#' @details
#' The function dispatches to internal plotting methods according to the
#' dimensionality of the fitted domain and the value of \code{mview}. For
#' two-dimensional triangulated domains, contour and surface plots are
#' supported. For three-dimensional tetrahedral domains, slice-based plots are
#' supported.
#'
#' Internally, the fitted model is evaluated on \code{Zgrid} through
#' prediction, so the resulting plot represents the fitted function on the
#' plotting grid rather than only at the observed training points.
#'
#' If a precomputed 3D array is supplied through \code{x$Yarray}, then
#' \code{plot.MPST(..., mview = "slice")} directly visualizes that array
#' without re-running prediction.
#'
#' @return
#' No value is returned. This function is called for its side effect of
#' generating a plot.
#'
#' @examples
#' \dontrun{
#' fit.obj <- fit.MPST(y ~ m(Z, V, Tr, d = 3, r = 1),
#'                     data = ex_square_train$m1$sigma01$tri50,
#'                     method = "G",
#'                     lambda = 10^seq(-6, 6, by = 0.5))
#'
#' plot(fit.obj, mview = "contour")
#' plot(fit.obj, mview = "surface")
#' plot(fit.obj, mview = "slice")
#' }
#'
#' @export
plot.MPST <- function(x, Zgrid = NULL, mview = NULL, slice_style = NULL, ...) {
  if (!inherits(x, "MPST")) {
    stop("Object must be of class 'MPST'")
  }
  
  if (is.null(mview)) {
    stop("mview must be specified (e.g., 'contour', 'surface', or 'slice').")
  }
  
  if (!mview %in% c("contour", "surface", "slice")) {
    stop("Invalid mview. Must be one of 'contour', 'surface', or 'slice'.")
  }
  
  has_precomputed_array <- !is.null(x$Yarray)
  
  if (has_precomputed_array) {
    if (mview != "slice") {
      stop("Precomputed Yarray objects currently support only mview = 'slice'.")
    }
    plot.slice.mpst(x, Zgrid = Zgrid, slice_style = slice_style, ...)
    return(invisible())
  }
  
  nd <- ncol(x$Tr)
  
  if ((mview == "contour") && (nd == 3)) {
    fig <- plot.contour.mpst(x, Zgrid = Zgrid, ...)
    print(fig)
  } else if ((mview == "surface") && (nd == 3)) {
    fig <- plot.surface.mpst(x, Zgrid = Zgrid, ...)
    print(fig)
  } else if ((mview == "slice") && (nd == 4)) {
    plot.slice.mpst(x, Zgrid = Zgrid, slice_style = slice_style, ...)
  } else {
    stop("Invalid mview or unsupported dimensionality")
  }
  
  invisible()
}


#' Initialize Grid for MPST Plotting
#'
#' @description
#' Generates a grid for plotting based on the input MPST model and dimensionality.
#' This function is used internally and is not intended for direct use by end users.
#'
#' @param mfit An MPST model fit object.
#' @param Zgrid An optional grid. If NULL, a grid is created automatically.
#' @param n1 Number of points in the first dimension.
#' @param n2 Number of points in the second dimension.
#' @param n3 Number of points in the third dimension (if applicable).
#'
#' @return A list containing the generated grid and corresponding coordinate vectors.
#'
#' @keywords internal
initialize.grid <- function(mfit, Zgrid = NULL, n1 = 101, n2 = 101, n3 = 101) {
  nd <- ncol(mfit$Tr)
  
  if (is.null(Zgrid)) {
    limits <- apply(mfit$Z[, 1:(nd - 1), drop = FALSE], 2, function(col) {
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
#' @description Generates a contour plot for a 2D MPST model. This function is used internally
#' by `plot.MPST` and is not intended for direct use by end users.
#'
#' @importFrom plotly plot_ly
#'
#' @param mfit An MPST model fit object.
#' @param Zgrid An optional grid for plotting. If NULL, a grid will be generated.
#'
#' @return A plotly object representing the contour plot.
#'
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
  
  z1 <- matrix(mpred$Ypred, nrow = length(u1), ncol = length(v1), byrow = FALSE)
  
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
#'
#' @importFrom plotly plot_ly add_trace layout
#'
#' @param mfit An MPST model fit object.
#' @param Zgrid An optional grid for plotting. If NULL, a grid is generated automatically.
#'
#' @return A plotly object representing the surface plot.
#'
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
  
  z1 <- matrix(mpred$Ypred, nrow = length(u1), ncol = length(v1), byrow = FALSE)
  df1 <- data.frame(x = mfit$Z[, 1], y = mfit$Z[, 2], z = mfit$Y)
  
  fig <- plotly::plot_ly(
    x = u1,
    y = v1,
    z = t(z1),
    type = "surface"
  )
  
  fig <- plotly::add_trace(
    fig,
    data = df1,
    x = ~x, y = ~y, z = ~z,
    mode = "markers", type = "scatter3d",
    marker = list(size = 1, color = "black")
  )
  
  fig <- plotly::layout(fig,
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
#' @description
#' Generates interactive 2D slices for a 3D MPST model. This function is used internally
#' by `plot.MPST` and is not intended for direct use by end users.
#'
#' It also supports direct slice plotting from a precomputed 3D array stored in
#' \code{mfit$Yarray}, together with \code{mfit$xgrid}, \code{mfit$ygrid}, and
#' \code{mfit$zgrid}.
#'
#' @importFrom manipulate manipulate slider
#' @importFrom fields image.plot
#'
#' @param mfit An MPST model fit object, or an MPST-like object containing
#' precomputed \code{Yarray}, \code{xgrid}, \code{ygrid}, and \code{zgrid}.
#' @param Zgrid An optional grid for plotting. If NULL, a grid is generated automatically.
#' @param slice_style A character string specifying the slice plotting style.
#' Supported values are \code{"default"} and \code{"brain"}.
#'
#' @return A manipulate object that provides interactive slices through the 3D array.
#'
#' @keywords internal
plot.slice.mpst <- function(mfit, Zgrid = NULL, slice_style = NULL) {
  if (!requireNamespace("manipulate", quietly = TRUE)) {
    stop("The 'manipulate' package is required for interactive plots.")
  }
  
  has_precomputed_array <- !is.null(mfit$Yarray)
  
  if (has_precomputed_array) {
    if (length(dim(mfit$Yarray)) != 3) {
      stop("mfit$Yarray must be a 3D array.")
    }
    if (is.null(mfit$xgrid) || is.null(mfit$ygrid) || is.null(mfit$zgrid)) {
      stop("Precomputed array objects must contain xgrid, ygrid, and zgrid.")
    }
    
    new.array <- mfit$Yarray
    z1.grid <- mfit$xgrid
    z2.grid <- mfit$ygrid
    z3.grid <- mfit$zgrid
    
    if (is.null(slice_style) && !is.null(mfit$slice_style)) {
      slice_style <- mfit$slice_style
    }
    if (is.null(slice_style)) {
      slice_style <- "default"
    }
    
  } else {
    if (!("Tr" %in% names(mfit)) || !("Z" %in% names(mfit))) {
      stop("The input object 'mfit' must contain 'Tr' and 'Z'.")
    }
    
    grid_info <- initialize.grid(mfit, Zgrid, n1 = 26, n2 = 28, n3 = 32)
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
    
    new.array <- array(
      mpred$Ypred,
      dim = c(length(z1.grid), length(z2.grid), length(z3.grid))
    )
    
    if (is.null(slice_style)) {
      slice_style <- "default"
    }
  }
  
  dim.size <- dim(new.array)
  
  select_color_palette <- function(color_choice, style = "default") {
    base_palette <- switch(as.character(color_choice),
                           "1" = gray.colors(64),
                           "2" = rainbow(64),
                           "3" = heat.colors(64),
                           "4" = terrain.colors(64),
                           "5" = topo.colors(64),
                           "6" = cm.colors(64)
    )
    
    if (identical(style, "brain")) {
      c("#FFFFFF00", base_palette)
    } else {
      base_palette
    }
  }
  
  plot_slices <- function(axial_slice, coronal_slice, sagittal_slice, color = 1) {
    par(mfrow = c(1, 3))
    
    if (identical(slice_style, "brain")) {
      arr_plot <- new.array
      arr_plot[is.na(arr_plot)] <- -1
      
      col_palette <- select_color_palette(color, style = "brain")
      zlim_use <- c(-1, max(arr_plot, na.rm = TRUE))
      
      fields::image.plot(
        1:dim.size[1], 1:dim.size[2], arr_plot[, , axial_slice],
        main = paste("Axial Plane (z =", axial_slice, ")"),
        xlab = "",
        ylab = "",
        col = col_palette,
        xlim = rev(c(1, dim.size[1])),
        ylim = c(-5, max(110, dim.size[2] + 5)),
        axes = FALSE,
        zlim = zlim_use,
        useRaster = TRUE
      )
      
      fields::image.plot(
        1:dim.size[1], 1:dim.size[3], arr_plot[, coronal_slice, ],
        main = paste("Coronal Plane (y =", coronal_slice, ")"),
        xlab = "",
        ylab = "",
        col = col_palette,
        xlim = rev(c(-5, max(80, dim.size[1]))),
        ylim = c(-60, max(130, dim.size[3] + 10)),
        axes = FALSE,
        zlim = zlim_use,
        useRaster = TRUE
      )
      
      fields::image.plot(
        1:dim.size[2], 1:dim.size[3], arr_plot[sagittal_slice, , ],
        main = paste("Sagittal Plane (x =", sagittal_slice, ")"),
        xlab = "",
        ylab = "",
        col = col_palette,
        xlim = rev(c(-5, max(100, dim.size[2]))),
        ylim = c(-70, max(150, dim.size[3] + 10)),
        axes = FALSE,
        zlim = zlim_use,
        useRaster = TRUE
      )
      
    } else {
      col_palette <- select_color_palette(color, style = "default")
      
      axial_mat <- new.array[, , axial_slice, drop = TRUE]
      coronal_mat <- new.array[, coronal_slice, , drop = TRUE]
      sagittal_mat <- new.array[sagittal_slice, , , drop = TRUE]
      
      if (all(is.na(axial_mat))) {
        warning("Axial slice contains only NA values.")
        plot.new()
        text(0.5, 0.5, "Axial slice contains only NA values.", cex = 1.5)
      } else {
        fields::image.plot(
          x = z1.grid,
          y = z2.grid,
          z = axial_mat,
          main = paste("Axial Plane (z =", round(z3.grid[axial_slice], 4), ")"),
          xlab = "x",
          ylab = "y",
          col = col_palette,
          axes = FALSE,
          useRaster = TRUE
        )
      }
      
      if (all(is.na(coronal_mat))) {
        warning("Coronal slice contains only NA values.")
        plot.new()
        text(0.5, 0.5, "Coronal slice contains only NA values.", cex = 1.5)
      } else {
        fields::image.plot(
          x = z1.grid,
          y = z3.grid,
          z = coronal_mat,
          main = paste("Coronal Plane (y =", round(z2.grid[coronal_slice], 4), ")"),
          xlab = "x",
          ylab = "z",
          col = col_palette,
          axes = FALSE,
          useRaster = TRUE
        )
      }
      
      if (all(is.na(sagittal_mat))) {
        warning("Sagittal slice contains only NA values.")
        plot.new()
        text(0.5, 0.5, "Sagittal slice contains only NA values.", cex = 1.5)
      } else {
        fields::image.plot(
          x = z2.grid,
          y = z3.grid,
          z = sagittal_mat,
          main = paste("Sagittal Plane (x =", round(z1.grid[sagittal_slice], 4), ")"),
          xlab = "y",
          ylab = "z",
          col = col_palette,
          axes = FALSE,
          useRaster = TRUE
        )
      }
    }
  }
  
  manipulate::manipulate(
    plot_slices(axial_slice, coronal_slice, sagittal_slice, color),
    axial_slice = manipulate::slider(
      1, max(1, dim.size[3]),
      initial = max(1, dim.size[3] %/% 2),
      label = "Axial Slice"
    ),
    coronal_slice = manipulate::slider(
      1, max(1, dim.size[2]),
      initial = max(1, dim.size[2] %/% 2),
      label = "Coronal Slice"
    ),
    sagittal_slice = manipulate::slider(
      1, max(1, dim.size[1]),
      initial = max(1, dim.size[1] %/% 2),
      label = "Sagittal Slice"
    ),
    color = manipulate::slider(
      1, 6,
      initial = 1,
      label = "Color Table (1-Gray, 2-Rainbow, 3-Heat, 4-Terrain, 5-Topo, 6-Cyan Magenta)"
    )
  )
  
  invisible()
}
