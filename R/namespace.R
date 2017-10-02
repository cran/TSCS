#' @import stats
#' @import ggplot2
#' @importFrom tseries adf.test
#' @importFrom rgl plot3d
#' @importFrom rgl title3d
#' @importFrom grDevices gray
NULL

#' @export
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


#' @export
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


#' @export
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


#' @export
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


#' @export
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


#' @export
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


#' @export
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


#' @export
### plotting - show spatial locations with or without missing observation in a given spatial domain ###
plot3D_NA <- function(newdata, xlab = NULL, ylab = NULL, zlab = NULL, title = NULL, cex = 3, color = "orange", colorNA = "blue"){
  
  ### obtain ID of spatial locations with missing observation ###
  NA_id <- which(is.na(newdata[,4])==TRUE);
  
  ### plotting - 'rgl' package is needed ###
  X <- newdata[,1];Y <- newdata[,2];Z <- newdata[,3];
  
  if (is.null(xlab)==TRUE) { xlab = names(newdata)[1]; }
  if (is.null(ylab)==TRUE) { ylab = names(newdata)[2]; }
  if (is.null(zlab)==TRUE) { zlab = names(newdata)[3]; }
  if (is.null(title)==TRUE) { title = paste("spatial(cross-section) data at time of",names(newdata)[4]); }
  
  rgl::plot3d(X[NA_id], Y[NA_id], Z[NA_id], col = colorNA, type = "p", size = cex, xlab = xlab, ylab = ylab, zlab = zlab);
  rgl::plot3d(X[-NA_id], Y[-NA_id], Z[-NA_id], col = color, type = "p", size = cex, add = TRUE);
  rgl::title3d(main = title);
}


#' @export
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


#' @export
plot3D_map <- function(newdata, xlab = NULL, ylab = NULL, zlab = NULL, title = NULL, cex = 9, colorNA = "white"){
  
  ### obtain ID of spatial locations with missing observation ###
  NA_id <- which(is.na(newdata[,4])==TRUE);
  
  ### plotting - 'rgl' package is needed ###
  X <- newdata[,1];Y <- newdata[,2];Z <- newdata[,3];
  indexT <- newdata[-NA_id,4];
  
  if (is.null(xlab)==TRUE) { xlab = names(newdata)[1]; }
  if (is.null(ylab)==TRUE) { ylab = names(newdata)[2]; }
  if (is.null(zlab)==TRUE) { zlab = names(newdata)[3]; }
  if (is.null(title)==TRUE) { title = paste("spatial(cross-section) data at time of",names(newdata)[4]); }
  
  color <- grDevices::gray((indexT - min(indexT))/(max(indexT) - min(indexT))); # mapped to [0,1]
  rgl::plot3d(X[NA_id], Y[NA_id], Z[NA_id], col = colorNA, type = "p", size = cex, xlab = xlab, ylab = ylab, zlab = zlab);
  rgl::plot3d(X[-NA_id], Y[-NA_id], Z[-NA_id], col = color, type = "p", size = cex, add = TRUE);
  rgl::title3d(main = title);
}


#' @export
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


#' @export
### compute appraisal indexes (RMSE and std) ###
appraisal_index <- function(est, true){
  N <- length(true);
  RMSE <- sqrt(sum((est - true)^2)/N); # root mean squared error (standard error)
  std <- sd(est - true); # standard deviation of error
  return(list(RMSE = RMSE, std = std))
}


