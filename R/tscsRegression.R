#' The First Step of TSCS for 2D Rectangular Grid System - Regression
#'
#' To implement TSCS spatial interpolation for a spatial domain that is a 2D rectangular grid system,
#'   the first step is obtaining regression coefficient matrix, which can be done
#'   by function \code{tscsRegression}. It is the prerequisite of TSCS interpolation process
#'   because the 'matrix' derived from historical spatio-temporal data is the initial value of
#'   the second step - estimating missing observations.
#'
#' @param data data frame; should contain these variables in order: X coordinate, Y coordinate and observations
#'   as time goes on. That is to say, each row should include X and Y coordinate first, and then a time series.
#'   This is the historical spatio-temporal data that you intend to analyze as the basis for
#'   interpolation later on in \code{tscsEstimate}.
#' @param h numeric; side length of the unit grid in X coordinate direction.
#' @param v numeric; side length of the unit grid in Y coordinate direction.
#' @param alpha numeric; specify the significance level for ADF test, to test if the time series of a group of
#'   spatial locations are cointegrated. (default: 0.05)
#'
#' @details
#' \itemize{
#'   \item The second step of TSCS spatial interpolation should be carried out by function \code{tscsEstimate},
#'   where you have to input the cross-section data or pure spatial data of a particular time point
#'   you have selected, with missing observations that you want to predict.
#'   \item For 3D rectangular grid system, the procedure of TSCS stays the same.
#'   Please see \code{tscsRegression3D} and \code{tscsEstimate3D}.
#'   \item Attentions:
#'   (1) Since TSCS is only capable of interpolation but not extrapolation, it is necessary to highlight the
#'   difference between interior spatial locations and system boundary. Function \code{plot_dif} can help.
#'   (2) NA value in historical spatio-temporal data \code{data} is not allowed. Please handle them beforehand
#'   (such as filling these NA values through spatio-temporal kriging).
#' }
#'
#' @return A list of 2 is returned, including:
#' \describe{
#'   \item{\code{coef_matrix}}{data frame; regression coefficient matrix to be used as input parameter of function
#'   \code{tscsEstimate} in the second step of TSCS interpolation.}
#'   \item{\code{percentage}}{numeric; percentage of cointegrated relationships, a measurement of the degree
#'   it satisfies the assumption of cointegrated system. It is highly affected by parameter \code{alpha},
#'   the significance level you have set. Explicitly, smaller \code{alpha} results in smaller \code{percentage}.}
#' }
#'
#' @seealso \code{\link{tscsEstimate}}, \code{\link{tscsRegression3D}}, \code{\link{plot_dif}}
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


##### obtain regression coefficients for all interior spatial locations #####
tscsRegression <- function(data, h, v, alpha = 0.05){

  ### select interior spatial locations ###
  XY_coord <- data[,1:2]; # coordinates of spatial locations
  X <- XY_coord[,1];Y <- XY_coord[,2];
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
  id <- which(index==1); # ID of interior spatial locations

  ### calculate regression coefficients - 'tseries' package is needed ###
  flag <- rep(NA,N); # for the use of recording if a group of spatial locations are cointegrated
  intercept = beta1 = beta2 = beta3 = beta4 = beta5 = beta6 = beta7 = beta8 = rep(NA,N);
  coef_matrix <- data.frame(X,Y,intercept,beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8);
  for (i in id) # NB: time consuming
  {
    Xp <- X[i];Yp <- Y[i];
    adp1 <- t(data[which((X==Xp)&(Y==Yp - v)),c(-1,-2)]); # transposition for regression - lm()
    adp2 <- t(data[which((X==Xp)&(Y==Yp + v)),c(-1,-2)]);
    adp3 <- t(data[which((X==Xp - h)&(Y==Yp)),c(-1,-2)]);
    adp4 <- t(data[which((X==Xp + h)&(Y==Yp)),c(-1,-2)]);
    adp5 <- t(data[which((X==Xp - h)&(Y==Yp - v)),c(-1,-2)]);
    adp6 <- t(data[which((X==Xp + h)&(Y==Yp + v)),c(-1,-2)]);
    adp7 <- t(data[which((X==Xp - h)&(Y==Yp + v)),c(-1,-2)]);
    adp8 <- t(data[which((X==Xp + h)&(Y==Yp - v)),c(-1,-2)]);
    reg <- lm(t(data[i,c(-1,-2)])~adp1 + adp2 + adp3 + adp4 + adp5 + adp6 + adp7 + adp8);
    error <- residuals(reg); # residuals
    adf.resid <- tseries::adf.test(error);
    flag[i] <- ifelse(adf.resid$p.value<=alpha,1,0); # set significance level alpha - 1 for cointegration and 0 for no
    coef_matrix$intercept[i] <- reg$coef[1];
    coef_matrix$beta1[i] <- reg$coef[2];
    coef_matrix$beta2[i] <- reg$coef[3];
    coef_matrix$beta3[i] <- reg$coef[4];
    coef_matrix$beta4[i] <- reg$coef[5];
    coef_matrix$beta5[i] <- reg$coef[6];
    coef_matrix$beta6[i] <- reg$coef[7];
    coef_matrix$beta7[i] <- reg$coef[8];
    coef_matrix$beta8[i] <- reg$coef[9];
  }
  percentage <- sum(flag[!is.na(flag)])/(N - sum(is.na(flag))); # percentage of cointegrated relationships
  coef_matrix <- coef_matrix[!is.na(coef_matrix[,3]),];
  return(list(coef_matrix = coef_matrix, percentage = percentage));
}


