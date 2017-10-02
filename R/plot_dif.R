#' Plot Interior Spatial Locations and System Boundary - 2D Map
#'
#' \code{plot_dif} differentiates boundary and interior spatial locations in a spatial domain (a collection of
#' spatial locations with their coordinates). Since TSCS method is only capable of interpolation but not
#' extrapolation, it is necessary to highlight the difference between interior spatial locations and system boundary.
#'
#' @param coords data frame; should only contain the two variables: X coordinate and Y coordinate. Each row uniquely
#'   denotes a spatial location. (coordinates must be numeric)
#' @param h numeric; side length of the unit grid in X coordinate direction.
#' @param v numeric; side length of the unit grid in Y coordinate direction.
#' @param xlab a label for the x axis, defaults to the name of X coordinate.
#' @param ylab a label for the y axis, defaults to the name of Y coordinate.
#' @param title a main title for the plot.
#' @param cex numeric; size of plotting point for each spatial location. (default: 1)
#'
#' @details
#'   \code{plot_dif} is exclusive to 2D rectangular grid system. Similarly, if you want to fathom how this package
#'   handles 3D rectangular grid system, please refer to \code{plot3D_dif}.
#'
#' @seealso \code{\link{plot3D_dif}}, \code{\link{plot_NA}}, \code{\link{plot_map}}
#'
#' @examples
#' \dontrun{
#' 
#' ## TSCS spatial interpolation procedure:
#' 
#' basis <- tscsRegression(data = data, h = 1, v = 1, alpha = 0.01); # regression
#' basis$percentage # see the percentage of cointegrated relationships
#' est <- tscsEstimate(matrix = basis$coef_matrix, newdata = newdata, h = 1, v = 1); # estimation
#' str(est)
#'
#' ## comparison of estimates and true values:
#' 
#' plot_compare(est = est$estimate[,3], true = true) # graphic comparison
#' index <- appraisal_index(est = est$estimate[,3], true = true); # RMSE & std
#' index
#'
#' ## data visualization:
#' 
#' plot_dif(data = data[,1:2], h = 1, v = 1) # differentiate boundary and interior spatial locations
#' plot_NA(newdata = newdata) # show spatial locations with missing value, for a cross-section data
#' plot_map(newdata = newdata) # plot the 2D spatial map, for a cross-section data
#' }



### plotting - differentiate boundary and interior spatial locations in a given spatial domain ###
plot_dif <- function(coords, h, v, xlab = NULL, ylab = NULL, title = NULL, cex = 1){

  ### select interior spatial locations ###
  X <- coords[,1];Y <- coords[,2];
  N <- length(X); # total number of spatial locations
  index <- rep(NA,N); # 0 for boundary point while 1 for interior point
  for (i in 1:N)
  {
    Xp <- X[i];Yp <- Y[i];
    j1 <- (length(which((X==Xp - h)&(Y==Yp)))!=0);
    j2 <- (length(which((X==Xp + h)&(Y==Yp)))!=0);
    j3 <- (length(which((X==Xp)&(Y==Yp - v)))!=0);
    j4 <- (length(which((X==Xp)&(Y==Yp + v)))!=0);
    j5 <- (length(which((X==Xp - h)&(Y==Yp - v)))!=0);
    j6 <- (length(which((X==Xp + h)&(Y==Yp + v)))!=0);
    j7 <- (length(which((X==Xp + h)&(Y==Yp - v)))!=0);
    j8 <- (length(which((X==Xp - h)&(Y==Yp + v)))!=0);
    index[i] <- as.numeric(j1&j2&j3&j4&j5&j6&j7&j8);
  }

  ### plotting - 'ggplot2' package is needed ###
  if (is.null(xlab)==TRUE) { xlab = names(coords)[1]; }
  if (is.null(ylab)==TRUE) { ylab = names(coords)[2]; }
  if (is.null(title)==TRUE) { title = "boundary and interior spatial locations"; }
  
  g <- as.factor(index + 1);
  ggplot(NULL, aes(x = X, y = Y, group = g, colour = g, shape = g)) + geom_point(cex = cex) +
    labs(x = xlab, y = ylab, title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(panel.background = element_rect(fill = "grey95", colour = "grey50")) +
    theme(legend.position = "bottom", legend.background = element_rect(fill = "grey95", colour = "grey50")) +
    theme(legend.key = element_rect(colour = "grey95")) +
    theme(legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10)) +
    scale_colour_discrete(name = "spatial locations", labels = c("boundary point","interior point")) +
    scale_shape_discrete(name = "spatial locations", labels = c("boundary point","interior point"))
}


