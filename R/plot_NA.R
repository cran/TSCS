#' Visualize the Spatial Distribution of Missing Observations - 2D Map
#'
#' \code{plot_NA} shows spatial locations with or without missing observation. It is plotted based on
#' the cross-section data of a given time point, which is also often extracted from spatio-temporal data.
#'
#' @param newdata data frame; should only contain the three variables in order: X coordinate, Y coordinate and observation.
#'   This is the cross-section data or pure spatial data of a particular time point you have selected,
#'   with missing observations that you want to predict. (coordinates must be numeric)
#' @param xlab a label for the x axis, defaults to the name of X coordinate.
#' @param ylab a label for the y axis, defaults to the name of Y coordinate.
#' @param title a main title for the plot.
#' @param cex numeric; size of plotting point for each spatial location. (default: 1)
#'
#' @details
#'   \code{plot_NA} is exclusive to 2D rectangular grid system. Similarly, if you want to fathom how this package
#'   handles 3D rectangular grid system, please refer to \code{plot3D_NA}.
#'
#' @seealso \code{\link{plot3D_NA}}, \code{\link{plot_map}}, \code{\link{plot_dif}}
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


### plotting - show spatial locations with or without missing observation in a given spatial domain ###
plot_NA <- function(newdata, xlab = NULL, ylab = NULL, title = NULL, cex = 1){

  ### obtain ID of spatial locations with missing observation ###
  NA_id <- which(is.na(newdata[,3])==TRUE);

  ### plotting - 'ggplot2' package is needed ###
  X <- newdata[,1];Y <- newdata[,2];
  
  if (is.null(xlab)==TRUE) { xlab = names(newdata)[1]; }
  if (is.null(ylab)==TRUE) { ylab = names(newdata)[2]; }
  if (is.null(title)==TRUE) { title = paste("spatial(cross-section) data at time of", names(newdata)[3]); }
  
  g <- as.factor(!is.na(newdata[,3]));
  ggplot(NULL, aes(x = X, y = Y, group = g, colour = g, shape = g)) + geom_point(cex = cex) +
    labs(x = xlab, y = ylab, title = title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    theme(panel.background = element_rect(fill = "grey95", colour = "grey50")) +
    theme(legend.position = "bottom", legend.background = element_rect(fill = "grey95", colour = "grey50")) +
    theme(legend.key = element_rect(colour = "grey95")) +
    theme(legend.title = element_text(face = "bold", size = 11), legend.text = element_text(size = 10)) +
    scale_colour_discrete(name = "observations", labels = c("missing","not missing")) +
    scale_shape_discrete(name = "observations", labels = c("missing","not missing"))
}


