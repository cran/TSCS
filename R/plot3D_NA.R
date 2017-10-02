#' Visualize the Spatial Distribution of Missing Observations - 3D Map
#'
#' \code{plot3D_NA} shows spatial locations with or without missing observation. It is plotted based on
#' the cross-section data of a given time point, which is also often extracted from spatio-temporal data.
#'
#' @param newdata data frame; should only contain the four variables in order: X coordinate, Y coordinate, Z coordinate
#'   and observation. This is the cross-section data or pure spatial data of a particular time point you have selected,
#'   with missing observations that you want to predict. (coordinates must be numeric)
#' @param xlab a label for the x axis, defaults to the name of X coordinate.
#' @param ylab a label for the y axis, defaults to the name of Y coordinate.
#' @param zlab a label for the z axis, defaults to the name of Z coordinate.
#' @param title a main title for the plot.
#' @param cex numeric; size of plotting point for each spatial location. (default: 3)
#' @param color colour to be used to fill the spatial locations. (default: "orange")
#' @param colorNA colour for denoting missing values/observations. (default: "blue")
#'
#' @details
#' \itemize{
#'   \item The resulting plot is interactive.
#'   \item \code{plot3D_NA} is exclusive to 3D rectangular grid system. Similarly, if you want to fathom how
#'   this package handles 2D rectangular grid system, please refer to \code{plot_NA}.
#' }
#'
#' @seealso \code{\link{plot_NA}}, \code{\link{plot3D_map}}, \code{\link{plot3D_dif}}
#'
#' @examples
#' \dontrun{
#' 
#' ## TSCS spatial interpolation procedure:
#'
#' basis <- tscsRegression3D(data = data, h1 = 3.75, h2 = 2.5, v = 5, alpha = 0.01);
#' basis$percentage
#' est <- tscsEstimate3D(matrix = basis$coef_matrix, newdata = newdata, h1 = 3.75, h2 = 2.5, v = 5);
#' str(est)
#'
#' ## comparison of estimates and true values:
#'
#' plot_compare(est = est$estimate[,4], true = true)
#' index <- appraisal_index(est = est$estimate[,4], true = true);
#' index
#'
#' ## data visualization:
#' 
#' plot3D_dif(data = data[,1:3], h1 = 3.75, h2 = 2.5, v = 5)
#' plot3D_NA(newdata = newdata)
#' plot3D_map(newdata = newdata)
#' }



### plotting - show spatial locations with or without missing observation in a given spatial domain ###
plot3D_NA <- function(newdata, xlab = NULL, ylab = NULL, zlab = NULL, title = NULL, cex = 3, color = "orange", colorNA = "blue"){

  ### obtain ID of spatial locations with missing observation ###
  NA_id <- which(is.na(newdata[,4])==TRUE);

  ### plotting - 'rgl' package is needed ###
  X <- newdata[,1];Y <- newdata[,2];Z <- newdata[,3];
  
  if (is.null(xlab)==TRUE) { xlab = names(newdata)[1]; }
  if (is.null(ylab)==TRUE) { ylab = names(newdata)[2]; }
  if (is.null(zlab)==TRUE) { zlab = names(newdata)[3]; }
  if (is.null(title)==TRUE) { title = paste("spatial(cross-section) data at time of",names(newdata)[4]); }
  
  rgl::plot3d(X[NA_id], Y[NA_id], Z[NA_id], col = colorNA, type = "p", size = cex, xlab = xlab, ylab = ylab, zlab = zlab);
  rgl::plot3d(X[-NA_id], Y[-NA_id], Z[-NA_id], col = color, type = "p", size = cex, add = TRUE);
  rgl::title3d(main = title);
}


