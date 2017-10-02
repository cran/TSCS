#' Plot Interior Spatial Locations and System Boundary - 3D Map
#'
#' \code{plot3D_dif} differentiates boundary and interior spatial locations in a spatial domain (a collection of
#' spatial locations with their coordinates). Since TSCS method is only capable of interpolation but not
#' extrapolation, it is necessary to highlight the difference between interior spatial locations and system boundary.
#'
#' @param coords data frame; should only contain the three variables: X coordinate, Y coordinate and Z coordinate.
#'   Each row uniquely denotes a spatial location. (coordinates must be numeric)
#' @param h1 numeric; side length of the unit cubic grid in X coordinate direction (horizontal).
#' @param h2 numeric; side length of the unit cubic grid in Y coordinate direction (horizontal).
#' @param v numeric; side length of the unit cubic grid in Z coordinate direction (vertical).
#' @param xlab a label for the x axis, defaults to the name of X coordinate.
#' @param ylab a label for the y axis, defaults to the name of Y coordinate.
#' @param zlab a label for the z axis, defaults to the name of Z coordinate.
#' @param title a main title for the plot.
#' @param cex numeric; size of point to be plotted for each spatial location. (default: 3)
#'
#' @details
#' \itemize{
#'   \item The resulting plot is interactive, where the red points are interior spatial locations
#'   while the black points denote system boundary.
#'   \item \code{plot3D_dif} is exclusive to 3D rectangular grid system. Similarly, if you want to fathom how
#'   this package handles 2D rectangular grid system, please refer to \code{plot_dif}.
#' }
#'
#' @seealso \code{\link{plot_dif}}, \code{\link{plot3D_NA}}, \code{\link{plot3D_map}}
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



### plotting - differentiate boundary and interior spatial locations in a given spatial domain ###
plot3D_dif <- function(coords, h1, h2, v, xlab = NULL, ylab = NULL, zlab = NULL, title = NULL, cex = 3){

  ### select interior spatial locations ###
  X <- coords[,1];Y <- coords[,2];Z <- coords[,3];
  N <- length(X); # total number of spatial locations
  index <- rep(NA,N); # 0 for boundary point while 1 for interior point
  for (i in 1:N)
  {
    Xp <- X[i];Yp <- Y[i];Zp <- Z[i];
    j1 <- (length(which((X==Xp - h1)&(Y==Yp)&(Z==Zp)))!=0);
    j2 <- (length(which((X==Xp + h1)&(Y==Yp)&(Z==Zp)))!=0);
    j3 <- (length(which((X==Xp)&(Y==Yp - h2)&(Z==Zp)))!=0);
    j4 <- (length(which((X==Xp)&(Y==Yp + h2)&(Z==Zp)))!=0);
    j5 <- (length(which((X==Xp)&(Y==Yp)&(Z==Zp - v)))!=0);
    j6 <- (length(which((X==Xp)&(Y==Yp)&(Z==Zp + v)))!=0);
    j7 <- (length(which((X==Xp + h1)&(Y==Yp + h2)&(Z==Zp + v)))!=0);
    j8 <- (length(which((X==Xp - h1)&(Y==Yp - h2)&(Z==Zp + v)))!=0);
    j9 <- (length(which((X==Xp + h1)&(Y==Yp - h2)&(Z==Zp + v)))!=0);
    j10 <- (length(which((X==Xp - h1)&(Y==Yp + h2)&(Z==Zp + v)))!=0);
    j11 <- (length(which((X==Xp + h1)&(Y==Yp + h2)&(Z==Zp - v)))!=0);
    j12 <- (length(which((X==Xp - h1)&(Y==Yp - h2)&(Z==Zp - v)))!=0);
    j13 <- (length(which((X==Xp + h1)&(Y==Yp - h2)&(Z==Zp - v)))!=0);
    j14 <- (length(which((X==Xp - h1)&(Y==Yp + h2)&(Z==Zp - v)))!=0);
    index[i] <- as.numeric(j1&j2&j3&j4&j5&j6&j7&j8&j9&j10&j11&j12&j13&j14);
  }

  ### plotting - 'rgl' package is needed ###
  if (is.null(xlab)==TRUE) { xlab = names(coords)[1]; }
  if (is.null(ylab)==TRUE) { ylab = names(coords)[2]; }
  if (is.null(zlab)==TRUE) { zlab = names(coords)[3]; }
  if (is.null(title)==TRUE) { title = "boundary and interior spatial locations"; }
  
  rgl::plot3d(X, Y, Z, col = index + 1, type = "p", size = cex, xlab = xlab, ylab = ylab, zlab = zlab);
  rgl::title3d(main = title);
}


