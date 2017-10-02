#' Visualize Spatial(Cross-Section) Data of a Given Time Point - 2D Map
#'
#' \code{plot_map} draws a two-dimensional spatial map. It is plotted based on the cross-section data
#'   of a given time point, which is also often extracted from spatio-temporal data.
#' 
#' @param newdata data frame; should only contain the three variables in order: X coordinate, Y coordinate and observation.
#'   This is the cross-section data or pure spatial data of a particular time point you have selected,
#'   with missing observations that you want to predict. (coordinates must be numeric)
#' @param xlab a label for the x axis, defaults to the name of X coordinate.
#' @param ylab a label for the y axis, defaults to the name of Y coordinate.
#' @param title a main title for the plot.
#' @param cex numeric; size of plotting point for each spatial locations. (default: 2)
#' @param shape either an integer specifying a symbol or a single character to be used as the default
#'   in plotting points. (default: 15)
#' @param low,high colours for low and high ends of the gradient. (default: "blue","red")
#' @param mid colour for midpoint of the gradient. (default: "yellow")
#' @param na.value colour for missing values/observations. (default: "white")
#' @param midpoint numeric; the midpoint of the gradient scale, defaults to the midpoint value of index presented.
#'
#' @details
#'   \code{plot_map} is exclusive to 2D rectangular grid system. Similarly, if you want to fathom how this package
#'   handles 3D rectangular grid system, please refer to \code{plot3D_map}.
#'
#' @seealso \code{\link{plot3D_map}}, \code{\link{plot_NA}}, \code{\link{plot_dif}}
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



plot_map <- function(newdata, xlab = NULL, ylab = NULL, title = NULL, cex = 2, shape = 15, low = "blue", mid = "yellow", high = "red", na.value = "white", midpoint = NULL){
  
  ### plotting - 'ggplot2' package is needed ###
  X <- newdata[,1];Y <- newdata[,2];
  index <- newdata[,3];
  
  if (is.null(xlab)==TRUE) { xlab = names(newdata)[1]; }
  if (is.null(ylab)==TRUE) { ylab = names(newdata)[2]; }
  if (is.null(title)==TRUE) { title = paste("spatial(cross-section) data at time of", names(newdata)[3]); }
  if (is.null(midpoint)==TRUE) { midpoint = 0.5*max(na.omit(index)) + 0.5*min(na.omit(index)); }
  
  ggplot(NULL, aes(x = X, y = Y)) + geom_point(cex = cex, shape = shape, aes(color = index)) +
    scale_colour_gradient2(low = low, mid = mid, high = high, na.value = na.value, midpoint = midpoint) +
    labs(x = xlab, y = ylab, title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(panel.background = element_rect(fill = "white", colour = "black")) +
    theme(legend.title = element_text(face = "bold", size = 11))
}


