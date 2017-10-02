#' The Second Step of TSCS for 3D Rectangular Grid System - Estimation
#'
#' \code{tscsEstimate} estimates the missing observations within the cross-section data (pure spatial data)
#'   of a particular time point you have selected, namely, the interpolation process.
#'
#' @param matrix data frame; the first return value \code{coef_matrix} of function \code{tscsRegression3D}
#'   in the first step of TSCS.
#' @param newdata data frame; should only contain the four variables in order: X coordinate, Y coordinate, Z coordinate
#'   and observation. This is the cross-section data or pure spatial data of a particular time point you have selected,
#'   with missing observations that you want to predict. (coordinates must be numeric)
#' @param h1 numeric; side length of the unit cubic grid in X coordinate direction (horizontal).
#' @param h2 numeric; side length of the unit cubic grid in Y coordinate direction (horizontal).
#' @param v numeric; side length of the unit cubic grid in Z coordinate direction (vertical).
#'
#' @details
#' \itemize{
#'   \item The first step of TSCS spatial interpolation should be carried out by function \code{tscsRegression3D},
#'   which is the prerequisite of \code{tscsEstimate3D}.
#'   \item For 2D rectangular grid system, the procedure of TSCS stays the same.
#'   Please see \code{tscsRegression} and \code{tscsEstimate}.
#'   \item Attentions:
#'   Since TSCS is only capable of interpolation but not extrapolation, please make sure that
#'   the missing observations in a given spatial domain are all located at interior spatial locations.
#'   Otherwise, extrapolation would occur with an error following.
#' }
#'
#' @return A list of 3 is returned, including:
#' \describe{
#'   \item{\code{estimate}}{data frame; estimate of missing observations which contains the 4 variables in order:
#'   X coordinate, Y coordinate, Z coordinate and estimation.}
#'   \item{\code{complete}}{data frame; an updated version of the cross-section data (pure spatial data) \code{newdata},
#'   with all of its missing observations interpolated.}
#'   \item{\code{NA_id}}{an integer vector; reveals the instance ID, in data frame \code{newdata},
#'   of spatial locations with missing observation.}
#' }
#'
#' @seealso \code{\link{tscsRegression3D}}, \code{\link{tscsEstimate}}, \code{\link{plot3D_NA}}, \code{\link{plot3D_map}}
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



##### estimate missing observations (interpolation) #####
tscsEstimate3D <- function(matrix, newdata, h1, h2, v){

  ### obtain ID of spatial locations with missing observation ###
  NA_id <- which(is.na(newdata[,4])==TRUE);

  ### estimate missing values - interpolation ###
  X <- newdata[,1];Y <- newdata[,2];Z <- newdata[,3];
  N <- length(NA_id); # number of missing values
  A <- matrix(0,N,N); # coefficients of the system of linear equations
  B <- rep(0,N);
  for (i in 1:N)
  {
    XN <- X[NA_id[i]];YN <- Y[NA_id[i]];ZN <- Z[NA_id[i]]; # X, Y and Z coordinates
    A[i,i] <- 1; # initial value
    B[i] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$intercept; # intercept - initial value
    if (!is.na(newdata[which((X==XN - h1)&(Y==YN)&(Z==ZN)),4])) {
      r1 <- newdata[which((X==XN - h1)&(Y==YN)&(Z==ZN)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta1*r1;    } else{
        s1 <- which(NA_id==which((X==XN - h1)&(Y==YN)&(Z==ZN)));
        A[i,s1] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta1*(-1); }
    if (!is.na(newdata[which((X==XN + h1)&(Y==YN)&(Z==ZN)),4])) {
      r2 <- newdata[which((X==XN + h1)&(Y==YN)&(Z==ZN)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta2*r2;    } else{
        s2 <- which(NA_id==which((X==XN + h1)&(Y==YN)&(Z==ZN)));
        A[i,s2] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta2*(-1); }
    if (!is.na(newdata[which((X==XN)&(Y==YN - h2)&(Z==ZN)),4])) {
      r3 <- newdata[which((X==XN)&(Y==YN - h2)&(Z==ZN)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta3*r3;    } else{
        s3 <- which(NA_id==which((X==XN)&(Y==YN - h2)&(Z==ZN)));
        A[i,s3] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta3*(-1); }
    if (!is.na(newdata[which((X==XN)&(Y==YN + h2)&(Z==ZN)),4])) {
      r4 <- newdata[which((X==XN)&(Y==YN + h2)&(Z==ZN)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta4*r4;    } else{
        s4 <- which(NA_id==which((X==XN)&(Y==YN + h2)&(Z==ZN)));
        A[i,s4] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta4*(-1); }
    if (!is.na(newdata[which((X==XN)&(Y==YN)&(Z==ZN - v)),4])) {
      r5 <- newdata[which((X==XN)&(Y==YN)&(Z==ZN - v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta5*r5;    } else{
        s5 <- which(NA_id==which((X==XN)&(Y==YN)&(Z==ZN - v)));
        A[i,s5] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta5*(-1); }
    if (!is.na(newdata[which((X==XN)&(Y==YN)&(Z==ZN + v)),4])) {
      r6 <- newdata[which((X==XN)&(Y==YN)&(Z==ZN + v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta6*r6;    } else{
        s6 <- which(NA_id==which((X==XN)&(Y==YN)&(Z==ZN + v)));
        A[i,s6] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta6*(-1); }
    if (!is.na(newdata[which((X==XN + h1)&(Y==YN + h2)&(Z==ZN + v)),4])) {
      r7 <- newdata[which((X==XN + h1)&(Y==YN + h2)&(Z==ZN + v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta7*r7;    } else{
        s7 <- which(NA_id==which((X==XN + h1)&(Y==YN + h2)&(Z==ZN + v)));
        A[i,s7] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta7*(-1); }
    if (!is.na(newdata[which((X==XN - h1)&(Y==YN - h2)&(Z==ZN + v)),4])) {
      r8 <- newdata[which((X==XN - h1)&(Y==YN - h2)&(Z==ZN + v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta8*r8;    } else{
        s8 <- which(NA_id==which((X==XN - h1)&(Y==YN - h2)&(Z==ZN + v)));
        A[i,s8] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta8*(-1); }
    if (!is.na(newdata[which((X==XN + h1)&(Y==YN - h2)&(Z==ZN + v)),4])) {
      r9 <- newdata[which((X==XN + h1)&(Y==YN - h2)&(Z==ZN + v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta9*r9;    } else{
        s9 <- which(NA_id==which((X==XN + h1)&(Y==YN - h2)&(Z==ZN + v)));
        A[i,s9] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta9*(-1); }
    if (!is.na(newdata[which((X==XN - h1)&(Y==YN + h2)&(Z==ZN + v)),4])) {
      r10 <- newdata[which((X==XN - h1)&(Y==YN + h2)&(Z==ZN + v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta10*r10;    } else{
        s10 <- which(NA_id==which((X==XN - h1)&(Y==YN + h2)&(Z==ZN + v)));
        A[i,s10] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta10*(-1); }
    if (!is.na(newdata[which((X==XN + h1)&(Y==YN + h2)&(Z==ZN - v)),4])) {
      r11 <- newdata[which((X==XN + h1)&(Y==YN + h2)&(Z==ZN - v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta11*r11;    } else{
        s11 <- which(NA_id==which((X==XN + h1)&(Y==YN + h2)&(Z==ZN - v)));
        A[i,s11] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta11*(-1); }
    if (!is.na(newdata[which((X==XN - h1)&(Y==YN - h2)&(Z==ZN - v)),4])) {
      r12 <- newdata[which((X==XN - h1)&(Y==YN - h2)&(Z==ZN - v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta12*r12;    } else{
        s12 <- which(NA_id==which((X==XN - h1)&(Y==YN - h2)&(Z==ZN - v)));
        A[i,s12] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta12*(-1); }
    if (!is.na(newdata[which((X==XN + h1)&(Y==YN - h2)&(Z==ZN - v)),4])) {
      r13 <- newdata[which((X==XN + h1)&(Y==YN - h2)&(Z==ZN - v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta13*r13;    } else{
        s13 <- which(NA_id==which((X==XN + h1)&(Y==YN - h2)&(Z==ZN - v)));
        A[i,s13] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta13*(-1); }
    if (!is.na(newdata[which((X==XN - h1)&(Y==YN + h2)&(Z==ZN - v)),4])) {
      r14 <- newdata[which((X==XN - h1)&(Y==YN + h2)&(Z==ZN - v)),4];
      B[i] <- B[i] + matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta14*r14;    } else{
        s14 <- which(NA_id==which((X==XN - h1)&(Y==YN + h2)&(Z==ZN - v)));
        A[i,s14] <- matrix[which((matrix$X==XN)&(matrix$Y==YN)&(matrix$Z==ZN)),]$beta14*(-1); }
  }
  solution <- solve(A,B); # estimate of missing observations
  estimate <- data.frame(X[NA_id], Y[NA_id], Z[NA_id], solution);
  names(estimate) <- c(names(newdata)[1], names(newdata)[2], names(newdata)[3], names(newdata)[4]);
  complete <- newdata;
  complete[NA_id,4] <- solution;
  names(complete) <- c(names(newdata)[1], names(newdata)[2], names(newdata)[3], names(newdata)[4]);
  return(list(estimate = estimate, complete = complete, NA_id = NA_id));
}


