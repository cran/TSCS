#' Graphic Comparison Between Estimates and True Values
#'
#' Provided that you have the true values of missing observations, you can compare them
#'   with the results of interpolation. \code{plot_compare} visualizes the comparison
#'   between estimates and true values. (NB: this plotting function can also be used
#'   in other similar situations involving comparison between estimates and true values.)
#'
#' @param est a numeric vector; estimations.
#' @param true a numeric vector; true values.
#' @param cex numeric; size of point to be plotted. (default: 1)
#' @param width numeric; width of fitted straight line. (default: 1)
#' @param P numeric, between 0 and 1; position for superimposing values of appraisal indexes. (default: 6/7)
#' @param AI logical; \code{TRUE} for presenting appraisal indexes while \code{FALSE} for not. (default: TRUE)
#'
#' @details
#' Attentions:
#' \itemize{
#'   \item The values in \code{est} and \code{true} vectors should be arranged in the same order,
#'   in correspondence with the sequence of observations.
#'   \item If the maximum value of either \code{est} or \code{true} is greater than 1000, or the
#'   minimum is smaller than -1000, please make appropriate transformation that limits your data
#'   to bound [-1000,1000].
#' }
#' 
#' In the plot:
#' \itemize{
#'   \item The big red point is the origin.
#'   \item The red line stands for straight line \code{y = x}.
#'   \item The blue line stands for fitted straight line.
#' }
#'
#' @seealso \code{\link{appraisal_index}}
#'
#' @examples
#' \dontrun{
#' 
#' ## TSCS spatial interpolation procedure:
#' 
#' basis <- tscsRegression(data = data, h = 1, v = 1, alpha = 0.01) # regression
#' basis$percentage # see the percentage of cointegrated relationships
#' est <- tscsEstimate(matrix = basis$coef_matrix, newdata = newdata, h = 1, v = 1) # estimation
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



### plotting - visualization of comparison ###
plot_compare <- function(est, true, cex = 1, width = 1, P = 6/7, AI = TRUE){

  print("The following warning message about 'Removed XXX rows ...' doesn't matter.");

  ### compute appraisal indexes ###
  N <- length(true);
  RMSE <- sqrt(sum((est - true)^2)/N); # RMSE
  std <- sd(est - true); # standard deviation of error

  ### format conversion ###
  RMSE <- as.character(round(RMSE,4));
  std <- as.character(round(std,4));
  joint <- paste("RMSE:", RMSE, "   ", "standard deviation of error:", std, seq = " ");

  ### set plotting range ###
  lowerX <- min(true) - (max(true) - min(true))/10;
  upperX <- max(true) + (max(true) - min(true))/10;
  lowerY <- min(est) - (max(est) - min(est))/10;
  upperY <- max(est) + (max(est) - min(est))/10;

  ### plotting - 'ggplot2' package is needed ###
  fit <- lm(est~true);
  magnitude <- floor(log10(0.5*mean(est) + 0.5*mean(true)));
  ss <- seq(-1200, 1200, 10^magnitude/20);
  
  if (AI==TRUE){
    ggplot(NULL, aes(x = true, y = est)) + geom_point(colour = "orange", cex = cex) +
      geom_line(aes(x = true, y = fitted(fit)), color = "blue", cex = width) +
      geom_point(aes(x = 0, y = 0), color = "red", cex = 3) +
      geom_line(aes(x = ss, y = ss), color = "red", cex = 0.5) + # maximum bound is [-1000,1000]
      xlim(lowerX,upperX) + ylim(lowerY,upperY) +
      labs(x = "true value", y = "estimate", title = "comparison between estimate and true value") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(panel.background = element_rect(fill = "grey95", colour = "grey50")) +
      annotate("text", x = (lowerX + upperX)*P, y = lowerY, label = joint, alpha = 0.8, size = 5) +
      annotate("text", x = lowerX, y = lowerY, label = "y = x", alpha = 0.8, size = 3.5, colour = "red3")
  } else{
    ggplot(NULL, aes(x = true, y = est)) + geom_point(colour = "orange", cex = cex) +
      geom_line(aes(x = true, y = fitted(fit)), color = "blue", cex = width) +
      geom_point(aes(x = 0, y = 0), color = "red", cex = 3) +
      geom_line(aes(x = ss, y = ss), color = "red", cex = 0.5) + # maximum bound is [-1000,1000]
      xlim(lowerX,upperX) + ylim(lowerY,upperY) +
      labs(x = "true value", y = "estimate", title = "comparison between estimate and true value") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      theme(panel.background = element_rect(fill = "grey95", colour = "grey50")) +
      annotate("text", x = lowerX, y = lowerY, label = "y = x", alpha = 0.8, size = 3.5, colour = "red3")
  }
}


