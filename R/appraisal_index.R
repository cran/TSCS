#' Compute Appraisal Index of Interpolation/Prediction Result
#'
#' Two appraisal indexes used for evaluating the result of interpolation/prediction - RMSE and
#'   standard deviation of error.
#'
#' @param est a numeric vector; estimations.
#' @param true a numeric vector; true values.
#'
#' @details
#' \itemize{
#'   \item The first appraisal index is RMSE, abbr. of root-mean-square error. It is used for measuring the differences
#'   between estimated values by a method and the values actually observed. Smaller RMSE means more accurate
#'   interpolation/prediction.
#'   \item The second appraisal index is standard deviation of error, which is used for measuring how far the errors
#'   are spread out from their mean, namely, stability of errors. Smaller value means greater stability of errors,
#'   suggesting that errors would not fluctuate heavily due to difference of data.
#' }
#'
#' @return A list of 2 is returned, including:
#' \describe{
#'   \item{\code{RMSE}}{numeric; RMSE.}
#'   \item{\code{std}}{numeric; standard deviation of error.}
#' }
#'
#' @seealso \code{\link{plot_compare}}
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



### compute appraisal indexes (RMSE and std) ###
appraisal_index <- function(est, true){
  N <- length(true);
  RMSE <- sqrt(sum((est - true)^2)/N); # root mean squared error (standard error)
  std <- sd(est - true); # standard deviation of error
  return(list(RMSE = RMSE, std = std))
}


