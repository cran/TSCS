#' The Second Step of TSCS for 2D Rectangular Grid System - Estimation
#'
#' \code{tscsEstimate} estimates the missing observations within the cross-section data (pure spatial data)
#'   of a particular time point you have selected, namely, the interpolation process.
#'
#' @param matrix data frame; the first return value \code{coef_matrix} of function \code{tscsRegression}
#'   in the first step of TSCS.
#' @param newdata data frame; should only contain the three variables in order: X coordinate, Y coordinate and observation.
#'   This is the cross-section data or pure spatial data of a particular time point you have selected,
#'   with missing observations that you want to predict. (coordinates must be numeric)
#' @param h numeric; side length of the unit grid in X coordinate direction.
#' @param v numeric; side length of the unit grid in Y coordinate direction.
#'
#' @details
#' \itemize{
#'   \item The first step of TSCS spatial interpolation should be carried out by function \code{tscsRegression},
#'   which is the prerequisite of \code{tscsEstimate}.
#'   \item For 3D rectangular grid system, the procedure of TSCS stays the same.
#'   Please see \code{tscsRegression3D} and \code{tscsEstimate3D}.
#'   \item Attentions:
#'   Since TSCS is only capable of interpolation but not extrapolation, please make sure that
#'   the missing observations in a given spatial domain are all located at interior spatial locations.
#'   Otherwise, extrapolation would occur with an error following.
#' }
#'
#' @return A list of 3 is returned, including:
#' \describe{
#'   \item{\code{estimate}}{data frame; estimate of missing observations which contains the 3 variables in order:
#'   X coordinate, Y coordinate and estimation.}
#'   \item{\code{complete}}{data frame; an updated version of the cross-section data (pure spatial data) \code{newdata},
#'   with all of its missing observations interpolated.}
#'   \item{\code{NA_id}}{an integer vector; reveals the instance ID, in data frame \code{newdata},
#'   of spatial locations with missing observation.}
#' }
#'
#' @seealso \code{\link{tscsRegression}}, \code{\link{tscsEstimate3D}}, \code{\link{plot_NA}}, \code{\link{plot_map}}
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



##### estimate missing observations (interpolation) #####
tscsEstimate <- function(matrix, newdata, h, v){

  ### obtain ID of spatial locations with missing observation ###
  NA_id <- which(is.na(newdata[,3])==TRUE);

  ### estimate missing values - interpolation ###
  X <- newdata[,1];Y <- newdata[,2];
  N <- length(NA_id); # number of missing values
  A <- matrix(0,N,N); # coefficients of the system of linear equations
  B <- rep(0,N);
  for (i in 1:N)
  {
    XN <- X[NA_id[i]];YN <- Y[NA_id[i]]; # X and Y coordinates
    A[i,i] <- 1; # initial value
    B[i] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$intercept; # intercept - initial value
    if (!is.na(newdata[which((X==XN)&(Y==YN - v)),3])) {
      r1 <- newdata[which((X==XN)&(Y==YN - v)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta1*r1;    } else{
        s1 <- which(NA_id==which((X==XN)&(Y==YN - v)));
        A[i,s1] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta1*(-1); }
    if (!is.na(newdata[which((X==XN)&(Y==YN + v)),3])) {
      r2 <- newdata[which((X==XN)&(Y==YN + v)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta2*r2;    } else{
        s2 <- which(NA_id==which((X==XN)&(Y==YN + v)));
        A[i,s2] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta2*(-1); }
    if (!is.na(newdata[which((X==XN - h)&(Y==YN)),3])) {
      r3 <- newdata[which((X==XN - h)&(Y==YN)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta3*r3;    } else{
        s3 <- which(NA_id==which((X==XN - h)&(Y==YN)));
        A[i,s3] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta3*(-1); }
    if (!is.na(newdata[which((X==XN + h)&(Y==YN)),3])) {
      r4 <- newdata[which((X==XN + h)&(Y==YN)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta4*r4;    } else{
        s4 <- which(NA_id==which((X==XN + h)&(Y==YN)));
        A[i,s4] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta4*(-1); }
    if (!is.na(newdata[which((X==XN - h)&(Y==YN - v)),3])) {
      r5 <- newdata[which((X==XN - h)&(Y==YN - v)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta5*r5;    } else{
        s5 <- which(NA_id==which((X==XN-h)&(Y==YN - v)));
        A[i,s5] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta5*(-1); }
    if (!is.na(newdata[which((X==XN + h)&(Y==YN + v)),3])) {
      r6 <- newdata[which((X==XN+h)&(Y==YN + v)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta6*r6;    } else{
        s6 <- which(NA_id==which((X==XN+h)&(Y==YN + v)));
        A[i,s6] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta6*(-1); }
    if (!is.na(newdata[which((X==XN - h)&(Y==YN + v)),3])) {
      r7 <- newdata[which((X==XN - h)&(Y==YN + v)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta7*r7;    } else{
        s7 <- which(NA_id==which((X==XN - h)&(Y==YN + v)));
        A[i,s7] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta7*(-1); }
    if (!is.na(newdata[which((X==XN + h)&(Y==YN - v)),3])) {
      r8 <- newdata[which((X==XN+h)&(Y==YN - v)),3];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta8*r8;    } else{
        s8 <- which(NA_id==which((X==XN + h)&(Y==YN - v)));
        A[i,s8] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)),]$beta8*(-1); }
  }
  solution <- solve(A,B); # estimate of missing observations
  estimate <- data.frame(X[NA_id], Y[NA_id], solution);
  names(estimate) <- c(names(newdata)[1], names(newdata)[2], names(newdata)[3]);
  complete <- newdata;
  complete[NA_id,3] <- solution;
  names(complete) <- c(names(newdata)[1], names(newdata)[2], names(newdata)[3]);
  return(list(estimate = estimate, complete = complete, NA_id = NA_id));
}


