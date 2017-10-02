#' The First Step of TSCS for 3D Rectangular Grid System - Regression
#'
#' To implement TSCS spatial interpolation for a spatial domain that is a 3D rectangular grid system,
#'   the first step is obtaining regression coefficient matrix, which can be done
#'   by function \code{tscsRegression3D}. It is the prerequisite of TSCS interpolation process
#'   because the 'matrix' derived from historical spatio-temporal data is the initial value of
#'   the second step - estimating missing observations.
#'
#' @param data data frame; should contain these variables in order: X coordinate, Y coordinate, Z coordinate and
#'   observations as time goes on. That is to say, each row should include X, Y and Z coordinate first, and then
#'   a time series. This is the historical spatio-temporal data that you intend to analyze as the basis for
#'   interpolation later on in \code{tscsEstimate3D}.
#' @param h1 numeric; side length of the unit cubic grid in X coordinate direction (horizontal).
#' @param h2 numeric; side length of the unit cubic grid in Y coordinate direction (horizontal).
#' @param v numeric; side length of the unit cubic grid in Z coordinate direction (vertical).
#' @param alpha numeric; specify the significance level for ADF test, to test if the time series of a group of
#'   spatial locations are cointegrated. (default: 0.05)
#'
#' @details
#' \itemize{
#'   \item The second step of TSCS spatial interpolation should be carried out by function \code{tscsEstimate3D},
#'   where you have to input the cross-section data or pure spatial data of a particular time point
#'   you have selected, with missing observations that you want to predict.
#'   \item For 2D rectangular grid system, the procedure of TSCS stays the same.
#'   Please see \code{tscsRegression} and \code{tscsEstimate}.
#'   \item Attentions:
#'   (1) Since TSCS is only capable of interpolation but not extrapolation, it is necessary to highlight the
#'   difference between interior spatial locations and system boundary. Function \code{plot3D_dif} can help.
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
#' @seealso \code{\link{tscsEstimate3D}}, \code{\link{tscsRegression}}, \code{\link{plot3D_dif}}
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



##### obtain regression coefficients for all interior spatial locations #####
tscsRegression3D <- function(data, h1, h2, v, alpha = 0.05){

  ### select interior spatial locations ###
  XYZ_coord <- data[,1:3]; # coordinates of spatial locations
  X <- XYZ_coord[,1];Y <- XYZ_coord[,2];Z <- XYZ_coord[,3];
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
  id <- which(index==1); # ID of interior spatial locations

  ### calculate regression coefficients ###
  flag <- rep(NA,N); # for the use of recording if a group of spatial locations are cointegrated
  intercept = beta1 = beta2 = beta3 = beta4 = beta5 = beta6 = rep(NA,N);
  beta7 = beta8 = beta9 = beta10 = beta11 = beta12 = beta13 = beta14 = rep(NA,N);
  coef_matrix <- data.frame(X,Y,Z,intercept,beta1,beta2,beta3,beta4,beta5,beta6);
  coef_matrix <- data.frame(coef_matrix,beta7,beta8,beta9,beta10,beta11,beta12,beta13,beta14);
  for (i in id) # NB: time consuming
  {
    Xp <- X[i];Yp <- Y[i];Zp <- Z[i];
    adp1 <- t(data[which((X==Xp - h1)&(Y==Yp)&(Z==Zp)),c(-1,-2,-3)]); # transposition for regression - lm()
    adp2 <- t(data[which((X==Xp + h1)&(Y==Yp)&(Z==Zp)),c(-1,-2,-3)]);
    adp3 <- t(data[which((X==Xp)&(Y==Yp - h2)&(Z==Zp)),c(-1,-2,-3)]);
    adp4 <- t(data[which((X==Xp)&(Y==Yp + h2)&(Z==Zp)),c(-1,-2,-3)]);
    adp5 <- t(data[which((X==Xp)&(Y==Yp)&(Z==Zp - v)),c(-1,-2,-3)]);
    adp6 <- t(data[which((X==Xp)&(Y==Yp)&(Z==Zp + v)),c(-1,-2,-3)]);
    adp7 <- t(data[which((X==Xp + h1)&(Y==Yp + h2)&(Z==Zp + v)),c(-1,-2,-3)]);
    adp8 <- t(data[which((X==Xp - h1)&(Y==Yp - h2)&(Z==Zp + v)),c(-1,-2,-3)]);
    adp9 <- t(data[which((X==Xp + h1)&(Y==Yp - h2)&(Z==Zp + v)),c(-1,-2,-3)]);
    adp10 <- t(data[which((X==Xp - h1)&(Y==Yp + h2)&(Z==Zp + v)),c(-1,-2,-3)]);
    adp11 <- t(data[which((X==Xp + h1)&(Y==Yp + h2)&(Z==Zp - v)),c(-1,-2,-3)]);
    adp12 <- t(data[which((X==Xp - h1)&(Y==Yp - h2)&(Z==Zp - v)),c(-1,-2,-3)]);
    adp13 <- t(data[which((X==Xp + h1)&(Y==Yp - h2)&(Z==Zp - v)),c(-1,-2,-3)]);
    adp14 <- t(data[which((X==Xp - h1)&(Y==Yp + h2)&(Z==Zp - v)),c(-1,-2,-3)]);
    reg <- lm(t(data[i,c(-1,-2,-3)])~adp1+adp2+adp3+adp4+adp5+adp6+adp7+adp8+adp9+adp10+adp11+adp12+adp13+adp14);
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
    coef_matrix$beta9[i] <- reg$coef[10];
    coef_matrix$beta10[i] <- reg$coef[11];
    coef_matrix$beta11[i] <- reg$coef[12];
    coef_matrix$beta12[i] <- reg$coef[13];
    coef_matrix$beta13[i] <- reg$coef[14];
    coef_matrix$beta14[i] <- reg$coef[15];
  }
  percentage <- sum(flag[!is.na(flag)])/(N - sum(is.na(flag))); # percentage of cointegrated relationships
  coef_matrix <- coef_matrix[!is.na(coef_matrix[,4]),];
  return(list(coef_matrix = coef_matrix, percentage = percentage));
}


