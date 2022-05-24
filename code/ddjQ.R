#####################################################################################################################
#
# Non-parametric Jump-Diffusion equation
# The code consists of two functions. 
# The first function, NWEst(), computes
# Kramers-Moyal (KM) coefficients using Nadaraya-Watson estimators for time series x. 
# The inputs and outputs of this function are as following.

# Inputs:
#  x0   : time series data
#  dx   : increment of x0 for smallest possible time step
#  nx   : number of dx to compute
#  dt   : time step
#  bw   : bandwidth for the kernel
#  na   : number of a values for computing the kernel
#  avec : mesh for the kernel

# Outputs:
#  K1     : 1st order KM coefficient or drift
#  K2     : 2nd order KM coefficient
#  D2     : diffusion
#  sigma2 : jump amplitude
#  lambda : jump rate
#  K4     : 4th order KM coefficient
#  K6     : 6th order KM coefficient

# The second function, ddj(), non-parametrically estimates  
# drift, diffusion and jump parameters of a time series.  
# It also interpolates to obtain the evolution of the nonparametric statistics in time. 
# The inputs and outputs of this function are as following.

# Inputs:
#  timeseries : numeric vector of the observed univariate timeseries values
#  bandwidth  : bandwidth of the kernel regressor (must be numeric). 
#  na         : number of points for computing the kernel (must be numeric)
#  dt         : time step

# Outputs:
#  avec   : mesh for which values of the nonparametric statistics are estimated
#  D1     : drift
#  D2     : diffusion
#  K2     : 2nd order KM coefficient
#  lambda : jump rate
#  sigma2 : jump amplitude
#  D2Jump : diffusion and jump ratio
#
#  Xvec1    : time series data as output
#  D1.t     : local drift
#  D2.t     : local diffusion
#  lambda.t : local jump rate
#  sigma2.t : local jump amplitude
#  D2Jump.t : local diffusion and jump ratio
#  Tvec1    : time index of time series data

######################################################################################################

NWEst <- function(x0,dx,nx,dt,bw,na,avec)  {
  if(!require(KernSmooth)){install.packages("KernSmooth")}
  require(KernSmooth)
  
  # Set up constants and useful preliminaries
  dt <- dt
  SF <- 1/(bw*sqrt(2*pi))  # scale factor for kernel calculation
  
  # Compute matrix of kernel values
  Kmat <- matrix(0,nrow=na,ncol=nx)
  for(i in 1:(nx)) {  # loop over columns (x0 values)
    Kmat[,i] <- SF*exp(-0.5*(x0[i]-avec)*(x0[i]-avec)/(bw*bw))
  }
  
  # Compute KM coefficients for each value of a
  K1 <- rep(0,na)
  K2 <- rep(0,na)
  K4 <- rep(0,na)
  K6 <- rep(0,na)
  
  for(i in 1:na) {  # loop over rows (a values)
    Ksum <- sum(Kmat[i,])  # sum of weights
    K1[i] <- (1/dt)*sum(Kmat[i,]*dx)/Ksum
    K2[i] <- (1/dt)*sum(Kmat[i,]*dx^2)/Ksum
    K4[i] <- (1/dt)*sum(Kmat[i,]*dx^4)/Ksum
    K6[i] <- (1/dt)*sum(Kmat[i,]*dx^6)/Ksum
  }
  
  # Compute jump parameters and diffusion
  sigma2 <- K6/(5*K4)
  lambda <- (K4/(3*sigma2^2))
  D2 <- (K2 -(lambda*sigma2))
  # set negative diffusion estimates to zero
  D2 <- ifelse(D2>0,D2,0)
  outlist <- list(K1 = K1, K2 = K2, D2 = D2, sigma2 = sigma2, 
                  lambda = lambda, K4 = K4, K6 = K6)
  return(outlist)
} 

#####################################################################################################################

ddj<-function(timeseries,bandwidth,na,dt,amin,amax){
  
  timeseries<-ts(timeseries) 
  
  d <- as.data.frame(timeseries)
  
  if (dim(timeseries)[2]==1){
    Y <- timeseries[,1]
    timeindex <- c(1:length(Y))*dt
  
  } else if (dim(timeseries)[2]==2){
    Y <- timeseries[,2]
    timeindex <- timeseries[,1]
  
  } else{
    print("not right format of timeseries input")
  }
  
  # Preliminaries
  Xvec1<-Y
  Tvec1<-timeindex
  dXvec1 <- diff(Y)
  dt <- dt
  #dt <- Tvec1[2]-Tvec1[1]
  bw <- bandwidth #*sd(Xvec1) # bandwidth 
  #amin <- min(Xvec1)
  amin <- amin
  #amax <- max(Xvec1)
  amax <- amax
  na <- na
  avec <- seq(amin,amax,length.out=na)
  nx <- length(dXvec1)
  
  ParEst <- NWEst(Xvec1,dXvec1,nx,dt,bw,na,avec)
  D1 <- ParEst$K1 #ParEst[[1]]
  K2 <- ParEst$K2 #ParEst[[2]]
  D2 <- ParEst$D2 #ParEst[[3]]
  sigma2 <- ParEst$sigma2 #ParEst[[4]]
  lambda <- ParEst$lambda #ParEst[[5]]
  K4 <- ParEst$K4 #ParEst[[6]]
  K6 <- ParEst$K6 #ParEst[[7]]

  D2Jump <- D2/(lambda*sigma2)
  
  D1.i <- approx(x=avec,y=D1,xout=Xvec1)
  D1.t <- D1.i$y 
  D2.i <- approx(x=avec,y=D2,xout=Xvec1)
  D2.t <- D2.i$y
  K4.i <- approx(x=avec,y=K4,xout=Xvec1)
  K4.t <- K4.i$y
  sigma2.i <- approx(x=avec,y=sigma2,xout=Xvec1)
  sigma2.t <- sigma2.i$y
  lambda.i <- approx(x=avec,y=lambda,xout=Xvec1)
  lambda.t <- lambda.i$y
  lambda.t <- ifelse(lambda.t < (1/dt), lambda.t, 1/dt)
  lambda <- ifelse(lambda < (1/dt), lambda, 1/dt)
  Jump.i <- approx(x=avec,y=lambda*sigma2,xout=Xvec1) 
  Jump.t <- Jump.i$y
  D2Jump.i <- approx(x=avec,y=D2Jump,xout=Xvec1)  
  D2Jump.t <- D2Jump.i$y

  Out <- list(avec = avec, D1 = D1, D2 = D2, K2 = K2, lambda = lambda, 
              sigma2 = sigma2, D2Jump = D2Jump, Xvec1 = Xvec1, 
	            D1.t = D1.t, D2.t = D2.t, K4.t = K4.t, sigma2.t = sigma2.t, 
	            lambda.t = lambda.t, Jump.t = Jump.t, D2Jump.t = D2Jump.t,
	            Tvec1 = Tvec1, K4 = K4, K6 = K6)
  return(Out)
}


#####################################################################################################################
# Calculate Q Criteria
# for different time lags

# Inputs:
#  x    : data to analyse for Q criteria
#  xbar : the value of x where Q should be evaluated
#  dx   : bin size around xbar
#  tau  : steps for time lag in sample

# Outputs:
#  Q     : Q value
#  eQ    : uncertainty of Q value
#####################################################################################################################

Q_maker <- function(x,xbar,dx,tau){

  Q <- eQ <- numeric(length = length(tau))
  
  for (i in 1:length(tau)){
    
    m4 <- m6 <- numeric(length = (length(x)-tau[i]))
    for (t in 1:(length(x)-tau[i])){
      if (x[t]>xbar-dx && x[t]<xbar+dx){
        m6[t] <- (x[t+tau[i]]-x[t])^6;
        m4[t] <- (x[t+tau[i]]-x[t])^4;
      }
    }
    Q[i]=mean(m6)/(5*mean(m4));
    eQ[i]=(sd(m6)+sd(m4))/sqrt(length(m4));
  }
  out <- list(Q = Q, eQ = eQ)
  return(out)
}


############################################################################################################################
#' Detrend a time series with a low-order polynomial.
#'
#' Remove low-order polynomial "trend" from a time series.  
#' @param x time series or vector of data values
#' @param n polynomial order
#' @return ts object with low-order LS polynomial removed.
#' @export
#' @examples
#' n <- 100
#' x <- arimaSim( n, ar = 0.9, int = 1 ) + 1:n / 20
#' plot( x )
#' d <- detrend( x, 1 )
#' plot( d )
#' 
#' detrend function from jrevenaugh/TSAUMN: Time Series Analysis Tools Used in ESCI5201 at the University of Minnesota
############################################################################################################################

detrend <- function( x, n = 1 ) {
  l <- length( x )
  r <- residuals( lm( x ~ stats::poly( 1:l, n ) ) )
  if ( is.ts( x ) ) return( ts( r, start = time( x )[1], deltat = deltat( x ) ) )
  return( r )
}
