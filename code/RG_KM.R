rm(list = ls())
graphics.off()

setwd("/home/pyeiphyolin/Documents/R_Code_Lin/snow/github/")
require(R.matlab)
require(plotrix)
require(KernSmooth)
require(gdata)
#install.packages("matrixStats")
library(matrixStats)
source("code/ddjQ.R")

d <- 5e-3          # diameter of the measuring tip
A <- 0.25*pi*d*d   # cross section of measuring tip
bw <- 0.3#1.7

# snow hardness data to analyze
data <- read.xls(xls ="data/Sample_profiles.xls")
z <- data$depth..mm.   # depth of the snow profile
#N <- length(z)  # length of the data
#PP1 <- data$X2017.01.09.Sample07..kPa. # in table as PP
#PP1o <- na.omit(PP1)
#PP2 <- data$X2017.01.31.Sample06..kPa.
#PP2o <- na.omit(PP2)
#DH1 <- data$X2017.01.31.Sample05..kPa. #
#DH1o <- na.omit(DH1)
#DH2 <- data$X2017.01.31.Sample11..kPa. #
#DH2o <- na.omit(DH2)
#DH3 <- data$X2017.01.31.Sample10..kPa.
#DH3o <- na.omit(DH3)
RG1 <- data$X2017.01.31.Sample01..kPa. #
RG1o <- na.omit(RG1)
RG2 <- data$X2017.01.31.Sample04..kPa. #
RG2o <- na.omit(RG2)
RG3 <- data$X2017.01.31.Sample03..kPa.
RG3o <- na.omit(RG3)
#RGlr1 <- data$X2017.01.09.Sample04..kPa. #
#RGlr1o <- na.omit(RGlr1)
#RGlr2 <- data$X2017.01.09.Sample05..kPa. #
#RGlr2o <- na.omit(RGlr2)
#RGlr3 <- data$X2017.01.31.Sample07..kPa. #
#RGlr3o <- na.omit(RGlr3)
#RGlr4 <- data$X2017.01.09.Sample01..kPa.
#RGlr4o <- na.omit(RGlr4)
#RGlr5 <- data$X2017.01.09.Sample02..kPa.
#RGlr5o <- na.omit(RGlr5)
#RGlr6 <- data$X2017.01.09.Sample03..kPa.
#RGlr6o <- na.omit(RGlr6)

######
# RG #
######

z1 <- z[1:length(RG1o)]
j <- which(z1>4)
zRG1 <- z1 <- z1[j]                  
RG1 <- RG1o[j]#/(1000*A)

z2 <- z[1:length(RG2o)]
j <- which(z2>4)
zRG2 <- z2 <- z2[j]                  
RG2 <- RG2o[j]#/(1000*A)

z3 <- z[1:length(RG3o)]
j <- which(z3>4)
zRG3 <- z3 <- z3[j]                  
RG3 <- RG3o[j]#/(1000*A)


# calculate the trend of profile using Gaussian kernel
RG1_tr <- ksmooth(z1, RG1, "normal", bandwidth = 0.6)$y   
RG2_tr <- ksmooth(z2, RG2, "normal", bandwidth = 0.6)$y
RG3_tr <- ksmooth(z3, RG3, "normal", bandwidth = 0.6)$y   

# detrend and normalize the snow profile
RG1 <- RG1 - RG1_tr
RG1 <- RG1/sd(RG1);
RG2 <- RG2 - RG2_tr
RG2 <- RG2/sd(RG2);
RG3 <- RG3 - RG3_tr
RG3 <- RG3/sd(RG3);

zRG_1 <- zRG1
RG_1 <- RG1
zRG_2 <- zRG2[1:750]
RG_2 <- RG2[1:750]
zRG_3 <- zRG2[751:1500]
RG_3 <- RG2[751:1500]
zRG_4 <- zRG3[1:750]
RG_4 <- RG3[1:750]
zRG_5 <- zRG3[751:1500]
RG_5 <- RG3[751:1500]

# store the detrended the snow profile in timeseries matrix
timeseries1 <- matrix(data = NA, nrow = length(zRG_1), ncol = 2)
timeseries1[,1] <- zRG_1
timeseries1[,2] <- RG_1

timeseries2 <- matrix(data = NA, nrow = length(zRG_2), ncol = 2)
timeseries2[,1] <- zRG_2
timeseries2[,2] <- RG_2

timeseries3 <- matrix(data = NA, nrow = length(zRG_3), ncol = 2)
timeseries3[,1] <- zRG_3
timeseries3[,2] <- RG_3

timeseries4 <- matrix(data = NA, nrow = length(zRG_4), ncol = 2)
timeseries4[,1] <- zRG_4
timeseries4[,2] <- RG_4

timeseries5 <- matrix(data = NA, nrow = length(zRG_5), ncol = 2)
timeseries5[,1] <- zRG_5
timeseries5[,2] <- RG_5

#bw <- 1.7#1.06*length(R)^(-0.2)#1.7  # bandwidth for Nadaraya-Watson estimator
na <- 40#20   # number of points for computing the kernel
dz <- z[2] - z[1]    # depth step
amin <- min(c(RG_1,RG_2,RG_3,RG_4,RG_5))
amax <- max(c(RG_1,RG_2,RG_3,RG_4,RG_5))

#source("ddj2.R")
# Estimate the drift, diffusion and jump parameters non-parametrically
ddjRG1 <- ddj(timeseries1, bandwidth=bw, na=na, dt = dz, amin, amax)
ddjRG2 <- ddj(timeseries2, bandwidth=bw, na=na, dt = dz, amin, amax)
ddjRG3 <- ddj(timeseries3, bandwidth=bw, na=na, dt = dz, amin, amax)
ddjRG4 <- ddj(timeseries4, bandwidth=bw, na=na, dt = dz, amin, amax)
ddjRG5 <- ddj(timeseries5, bandwidth=bw, na=na, dt = dz, amin, amax)

Rbin_RG <- Rbin <- ddjRG1$avec

D1mat <- cbind(ddjRG1$D1, ddjRG2$D1, ddjRG3$D1, ddjRG4$D1, ddjRG5$D1)
D1_RG <- D1 <- rowMeans(D1mat)
eD1_RG <- eD1 <- rowSds(D1mat)/sqrt(6)

K2mat <- cbind(ddjRG1$K2, ddjRG2$K2, ddjRG3$K2, ddjRG4$K2, ddjRG5$K2)
K2_RG <- K2 <- rowMeans(K2mat)
eK2_RG <- eK2 <- rowSds(K2mat)/sqrt(6)

K4mat <- cbind(ddjRG1$K4, ddjRG2$K4, ddjRG3$K4, ddjRG4$K4, ddjRG5$K4)
K4_RG <- K4 <- rowMeans(K4mat)
eK4_RG <- eK4 <- rowSds(K4mat)/sqrt(6)

K6mat <- cbind(ddjRG1$K6, ddjRG2$K6, ddjRG3$K6, ddjRG4$K6, ddjRG5$K6)
K6_RG <- K6 <- rowMeans(K6mat)
eK6_RG <- eK6 <- rowSds(K6mat)/sqrt(6)

sigma2_mat <- K6mat/(5*K4mat)
sigma2_RG <- sigma2 <- K6/(5*K4)
#esigma2_RG <- esigma2 <- sigma2*sqrt((eK6/K6)^2 + (eK4/K4)^2)
esigma2_RG <- esigma2 <- rowSds(sigma2_mat)/sqrt(6)

lambda_mat <- K4mat/(3*sigma2_mat*sigma2_mat)
lambda <- K4/(3*sigma2*sigma2)
#elambda <- lambda*sqrt((3*eK4/K4)^2 + (2*eK6/K6)^2)
elambda <- rowSds(lambda_mat)/sqrt(6)
lambda_RG <- lambda <- ifelse(lambda < (1/dz), lambda, 1/dz)
elambda_RG <- elambda <- ifelse(elambda < (1/dz), elambda, 1/dz)

D2mat <- K2mat - lambda_mat*sigma2_mat
D2 <- K2 - lambda*sigma2
#eD2_RG <- eD2 <- sqrt(eK2^2 + (eK4*(10/3)*(K4/K6))^2 + (eK6*(5/3)*(K4/K6)^2)^2)
eD2_RG <- eD2 <- rowSds(D2mat)/sqrt(6)
D2_RG <- D2 <- ifelse(D2>0,D2,0)

DiffJump_mat <- D2/(lambda_mat*sigma2_mat)
DiffJump_RG <- DiffJump <- D2/(lambda*sigma2)
#eDiffJump_RG <- eDiffJump <- sqrt((eK2*(3/5)*(K6/K4^2))^2 + (eK4*(6/5)*(K2*K6/K4^3))^2 + (eK6*(3/5)*(K2/K4^2))^2)
eDiffJump_RG <- eDiffJump <- rowSds(DiffJump_mat)/sqrt(6)

gamma = as.numeric(coef(lm(D1[Rbin<2 & Rbin>-2]~Rbin[Rbin<2 & Rbin>-2] + 0)))
1/gamma
lambda_mean = mean(lambda[Rbin<2 & Rbin>-2])
1/lambda_mean
lambda_0 = lambda[Rbin == min(abs(Rbin))]
1/lambda_0
sigma2_mean = mean(sigma2[Rbin<2 & Rbin>-2])
sigma2_mean
DiffJump_mean = mean(DiffJump[Rbin<2 & Rbin>-2])
DiffJump_mean

pdf('plots/D1_RG_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(Rbin,D1,type='p',lwd=1,col='black',
     xlab=expression(italic(R)^"'"),ylab=expression(italic(D)^"(1)"), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-4,4), ylim = c(-600,600))
abline(h = 0, lty=3 , col = 'black')
polygon(c(rev(Rbin), Rbin), c(rev(-eD1+D1), (eD1+D1)), col = 'grey80', border = NA)
points(Rbin,D1)
dev.off()

pdf('plots/D2_RG_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(Rbin,D2,type='p',lwd=1,col='black',
     xlab=expression(italic(R)^"'"),ylab=expression(italic(D)^"(2)"), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-4,4), ylim = c(0,200))
polygon(c(rev(Rbin), Rbin), c(rev(-eD2+D2), (eD2+D2)), col = 'grey80', border = NA)
points(Rbin,D2)
dev.off()

pdf('plots/K4_RG_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(Rbin,K4,type='p',lwd=1,col='black',
     xlab=expression(italic(R)^"'"),ylab=expression(italic(K)^"(4)"), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-4,4), ylim = c(0,1000))
polygon(c(rev(Rbin), Rbin), c(rev(-eK4+K4), (eK4+K4)), col = 'grey80', border = NA)
points(Rbin,K4)
dev.off()

pdf('plots/JumpProb_RG_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(Rbin,lambda*dz,type='p',lwd=1,col='black',
     xlab=expression(italic(R)^"'"),ylab=expression(italic(lambda~Delta*z)), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-4,4), ylim = c(0,1.1))
polygon(c(rev(Rbin), Rbin), c(rev(-elambda+lambda)*dz, (elambda+lambda)*dz), col = 'grey80', border = NA)
points(Rbin,lambda*dz)
abline(h = lambda_mean*dz, col = 'blue')
dev.off()

pdf('plots/JumpAmp_RG_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(Rbin,sigma2,type='p',lwd=1,col='black',
     xlab=expression(italic(R)^"'"),ylab=expression(italic(sigma)[italic(xi)]^2), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-4,4), ylim = c(0,3))
polygon(c(rev(Rbin), Rbin), c(rev(-esigma2+sigma2), (esigma2+sigma2)), col = 'grey80', border = NA)
points(Rbin,sigma2)
abline(h = sigma2_mean, col = 'blue')
dev.off()

pdf('plots/DiffJump_RG_.pdf',width = 4, height = 3.5)
par(mfrow=c(1,1),mar=c(4, 5, 0, 0),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(Rbin,DiffJump,type='p',lwd=1,col='black',
     xlab=expression(italic(R)^"'"),ylab=expression(italic(D)^"(2)"/italic(lambda~sigma)[italic(xi)]^2), 
     cex.lab = 1.3, cex.axis=1.3, xlim = c(-4,4), ylim = c(0,1.5)) #max(D2sig2_RG)))
polygon(c(rev(Rbin), Rbin), c(rev(-eDiffJump+DiffJump), (eDiffJump+DiffJump)), col = 'grey80', border = NA)
points(Rbin,DiffJump)
abline(h = DiffJump_mean, col = 'blue')
dev.off()

