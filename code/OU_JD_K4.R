rm(list=ls())
graphics.off()

library(plotrix)
library(tictoc)
tic()

setwd("/home/pyeiphyolin/Documents/R_Code_Lin/snow/github/")

# Number of data points and initial conditions
N <- 1e6#2e6
x0 <- 0
X <- x0
Xvec <- rep(NA,N) 
Xvec[1] <- x0
XD <- x0
XDvec <- rep(NA,N) 
XDvec[1] <- x0
XJ <- x0
XJvec <- rep(NA,N) 
XJvec[1] <- x0


# parameters for OU process with jump
# dX = -gamma*dt + sqrt(D)*dWt + xi*dJt

gamma <- 100#10
D <- 10#0.1
sigma2 <- 1#10#100
lambda <- 100#60#400#200#10 
dt <- 1e-3
#dt <- seq(from=0, to =1e-2, by = 4e-4)
#dt <- dt[2:length(dt)]

#Q <- QD <-  rep(NA,length(dt)) 
ldt <- lambda*dt

set.seed(1)
eta <- rnorm(N,0,1)
set.seed(2)
xi <- rnorm(N,0,sqrt(sigma2))
set.seed(3)
u <- runif(N)

ul <- (1-ldt)/2
ur <- (1+ldt)/2

Jvec <- rep(NA,N) 
Jvec[1] <- 0
no = 1

for (i in c(2:N)){
  
  if (u[i] >= ul & u[i] <= ur){
    J <- 1
  }else{
    J <- 0
  }
  Jvec[i] <- J
  XD = XD - (gamma * XD * dt) + (sqrt(dt*D)*eta[i])
  X =  X - (gamma * X * dt) + (sqrt(dt*D)*eta[i]) + (xi[i] * J)
  XJ =  XJ - (gamma * XJ * dt) + (xi[i] * J)
  
  XDvec[i] <- XD
  Xvec[i] <- X
  XJvec[i] <- XJ
  
}

t <- c(1:N)*dt
#Xvec_norm <- Xvec/sd(Xvec)
toc()

Xvec <- Xvec - mean(Xvec)
Xvecn <- Xvec/sd(Xvec)
XJvec <- XJvec - mean(XJvec)
XJvecn <- XJvec/sd(XJvec)
XDvec <- XDvec - mean(XDvec)
XDvecn <- XDvec/sd(XDvec)

ii <- c(1:400)
plot(t[ii], XDvec[ii], type = 'l', ylim = c(-5,5))
plot(t[ii], XJvec[ii], type = 'l', ylim = c(-5,5))
plot(t[ii], Xvec[ii], type = 'l', ylim = c(-5,5))

source("/home/pyeiphyolin/Documents/R_Code_Lin/snow/github/code/ddjQ.R")

XX <- data.frame(t, XDvec, XJvec, Xvec)
write.table(XX,"data/OU_J_D.txt")

#x <- Xvec
bw <- 0.3
na <- 51
avec <- seq(min(Xvecn), max(Xvecn), length.out = na)
avecJ <- seq(min(XJvecn), max(XJvecn), length.out = na)
avecD <- seq(min(XDvecn), max(XDvecn), length.out = na)

dX <- diff(Xvecn)
dXJ <- diff(XJvecn)
dXD <- diff(XDvecn)
nx <- length(dX)

library(tictoc)
tic()
KM <- NWEst(Xvecn,dX,nx,dt,bw,na,avec)
KMJ <- NWEst(XJvecn,dXJ,nx,dt,bw,na,avecJ)
KMD <- NWEst(XDvecn,dXD,nx,dt,bw,na,avecD)
toc()

K2 <- KM$K2
K4 <- KM$K4
lambda <- KM$lambda
sigma2 <- KM$sigma2
ltau <- lambda*dt
D2 <- K2 - lambda*sigma2
D2 <- ifelse(D2<0,0,D2)
diffJump <- D2/(lambda*sigma2)

K2J <- KMJ$K2
K4J <- KMJ$K4
lambdaJ <- KMJ$lambda
sigma2J <- KMJ$sigma2
ltauJ <- lambdaJ*dt
D2J <- K2J - lambdaJ*sigma2J
D2J <- ifelse(D2J<0,0,D2J)
diffJumpJ <- D2J/(lambdaJ*sigma2J)

K2D <- KMD$K2
K4D <- KMD$K4
lambdaD <- KMD$lambda
sigma2D <- KMD$sigma2
ltauD <- lambdaD*dt
D2D <- K2D - lambdaD*sigma2D
D2D <- ifelse(D2D<0,0,D2D)
diffJumpD <- D2D/(lambdaD*sigma2D)

plot(avec,K4, xlim = c(-4,4), ylim = c(0,3000))
points(avecD,K4D,pch=2)
points(avecJ,K4J,pch=3)

plot(avec,K2, xlim = c(-2,2), ylim = c(0,400))
points(avecD,K2D,pch=2)
points(avecJ,K2J,pch=3)

plot(avec,diffJump, xlim = c(-2,2), ylim = c(0,2))
points(avecD,diffJumpD,pch=2)
points(avecJ,diffJumpJ,pch=3)

ii <- c(1:400)
pdf('plots/OU_D_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(t[ii], XDvec[ii], type = 'l',
     xlab=expression(italic(t)~"[s]"),ylab=expression(italic(x)), 
     ylim = c(-4,4),
     cex.lab = 1.3, cex.axis=1.3)
dev.off()

pdf('plots/OU_J_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(t[ii], XJvec[ii], type = 'l',
     xlab=expression(italic(t)~"[s]"),ylab=expression(italic(x)), 
     ylim = c(-4,4),
     cex.lab = 1.3, cex.axis=1.3)
dev.off()

pdf('plots/OU_JD_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(t[ii], Xvec[ii], type = 'l',
     xlab=expression(italic(t)~"[s]"),ylab=expression(italic(x)), 
     ylim = c(-4,4),
     cex.lab = 1.3, cex.axis=1.3)
dev.off()


pdf('plots/K4_OU_D.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(avecD,K4D,type='p',lwd=1,col='black',
     xlab=expression(italic(x)),ylab=expression(italic(K)^"(4)"), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-2,2), ylim = c(0,2000))
dev.off()

pdf('plots/K4_OU_J.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(avecJ,K4J,type='p',lwd=1,col='black',
     xlab=expression(italic(x)),ylab=expression(italic(K)^"(4)"), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-2,2), ylim = c(0,2000))
dev.off()

pdf('plots/K4_OU_JD.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(avec,K4,type='p',lwd=1,col='black',
     xlab=expression(italic(x)),ylab=expression(italic(K)^"(4)"), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-2,2), ylim = c(0,2000))
dev.off()



#########################################################################

rm(list=ls())
graphics.off()

library(plotrix)
library(tictoc)
tic()

setwd("/home/pyeiphyolin/Documents/R_Code_Lin/snow/github/")

# Number of data points and initial conditions
N <- 1e6#2e6
x0 <- 0
X <- x0
Xvec <- rep(NA,N) 
Xvec[1] <- x0
XD <- x0
XDvec <- rep(NA,N) 
XDvec[1] <- x0
XJ <- x0
XJvec <- rep(NA,N) 
XJvec[1] <- x0


# parameters for OU process with jump
# dX = -gamma*dt + sqrt(D)*dWt + xi*dJt

gamma <- 100#10
D <- 5#10#5#0.1
sigma2 <- 1#0.5#1#10#100
lambda <- 100#60#400#200#10 
dt <- 1e-3
#dt <- seq(from=0, to =1e-2, by = 4e-4)
#dt <- dt[2:length(dt)]

#Q <- QD <-  rep(NA,length(dt)) 
ldt <- lambda*dt

set.seed(1)
eta <- rnorm(N,0,1)
set.seed(2)
xi <- rnorm(N,0,sqrt(sigma2))
set.seed(3)
u <- runif(N)

ul <- (1-ldt)/2
ur <- (1+ldt)/2

Jvec <- rep(NA,N) 
Jvec[1] <- 0
no = 1

for (i in c(2:N)){
  
  if (u[i] >= ul & u[i] <= ur){
    J <- 1
  }else{
    J <- 0
  }
  Jvec[i] <- J
  XD = XD - (gamma * XD * dt) + (sqrt(dt*D)*eta[i])
  X =  X - (gamma * X * dt) + (sqrt(dt*D)*eta[i]) + (xi[i] * J)
  XJ =  XJ - (gamma * XJ * dt) + (xi[i] * J)
  
  XDvec[i] <- XD
  Xvec[i] <- X
  XJvec[i] <- XJ
  
}

t <- c(1:N)*dt
#Xvec_norm <- Xvec/sd(Xvec)
toc()

Xvec <- Xvec - mean(Xvec)
Xvec <- Xvec/sd(Xvec)
XJvec <- XJvec - mean(XJvec)
XJvec <- XJvec/sd(XJvec)
XDvec <- XDvec - mean(XDvec)
XDvec <- XDvec/sd(XDvec)

ii <- c(1:400)
#plot(t[ii], XDvec[ii], type = 'l', ylim = c(-5,5))
#plot(t[ii], XJvec[ii], type = 'l', ylim = c(-5,5))
plot(t[ii], Xvec[ii], type = 'l', ylim = c(-5,5))

#x <- Xvec
bw <- 0.3
na <- 51
avec <- seq(min(Xvec), max(Xvec), length.out = na)
avecJ <- seq(min(XJvec), max(XJvec), length.out = na)
avecD <- seq(min(XDvec), max(XDvec), length.out = na)

dX <- diff(Xvec)
dXJ <- diff(XJvec)
dXD <- diff(XDvec)
nx <- length(dX)

source("/home/pyeiphyolin/Documents/R_Code_Lin/snow/github/code/ddjQ.R")
library(tictoc)
tic()
KM <- NWEst(Xvec,dX,nx,dt,bw,na,avec)
toc()

K2 <- KM$K2
K4 <- KM$K4
lambda <- KM$lambda
sigma2 <- KM$sigma2
ltau <- lambda*dt
D2 <- K2 - lambda*sigma2
D2 <- ifelse(D2<0,0,D2)
diffJump <- D2/(lambda*sigma2)

ii <- c(1:400)
pdf('plots/OU_0_05.pdf',width = 4, height = 3.5)
#pdf('plots/OU_0_1.pdf',width = 4, height = 3.5)
#pdf('plots/OU_0_2.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(t[ii], Xvec[ii], type = 'l',
     xlab=expression(italic(t)~"[s]"),ylab=expression(italic(x)), 
     ylim = c(-4,4),
     cex.lab = 1.3, cex.axis=1.3)
dev.off()

pdf('plots/diffjump_OU_0_05.pdf',width = 4, height = 3.5)
#pdf('plots/diffjump_OU_0_1.pdf',width = 4, height = 3.5)
#pdf('plots/diffjump_OU_0_2.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(avec,diffJump,type='p',lwd=1,col='black',
     xlab=expression(italic(x)),ylab=expression(italic(D)^"(2)"/italic(lambda~sigma)[italic(xi)]^2), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-2,2), ylim = c(0,0.5))
abline(h = 0.05, col = 'blue')
#abline(h = 0.1, col = 'blue')
#abline(h = 0.2, col = 'blue')
dev.off()
