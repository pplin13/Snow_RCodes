rm(list = ls())
graphics.off()
setwd("/home/pyeiphyolin/Documents/R_Code_Lin/snow/github/")
require(R.matlab)
require(KernSmooth)
source('code/ddjQ.R')
data <- read.csv("data/SMP_profiles.csv")
z = data$depth..mm.
R = data$S29M0275..kPa.#data$S29M0277..kPa.
#ind = which(z >= 600)#which(z >= 600 & z < 800)
ind = c(1:300001)
z = z[ind]
R = R[ind]
#rm(data)

plot(z,R,type = 'l')

#R_tr = ksmooth(z,R, "normal", bandwidth = 0.6)$y
#R_detr = R - R_tr
#R_detr = R_detr/sd(R_detr)

bw = 0.3#1.7
na=40
dz = z[2] - z[1]
dw = 500#1000;
Lwindow = 500#1000;
Nwindow = length(seq(from = 1, to = (length(R) - Lwindow + 1), by = dw));

z_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
R_detr_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
R_tr_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
R_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
Qrt_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
D1t_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
K4t_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
D2Jumpt_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
lambdat_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
Jumpt_window = matrix(NA, nrow = Nwindow, ncol = Lwindow)
i = 0;

for (j in seq(from = 1, to = (length(R) - Lwindow + 1), by = dw)){
  
  i = i+1
  R_tr = ksmooth(z[j:(j-1+Lwindow)],R[j:(j-1+Lwindow)], "normal", bandwidth = 0.6)$y
  R_detr = R[j:(j-1+Lwindow)] - R_tr
  R_detr = R_detr/sd(R_detr)
  amin = min(R_detr)
  amax = max(R_detr)
  z_window[i,] = z[j:(j-1+Lwindow)]
  R_detr_window[i,] = R_detr
  R_tr_window[i,] = R_tr
  Rts <- cbind(z[j:(j-1+Lwindow)],R_detr)
  #Qr_window[i,] <- Qrt(Rts,windows=5,bandwidth=0.3,na=na,DT=dz,logtransform=FALSE,interpolate=FALSE)
  ddj.t <- ddj(Rts,bandwidth = bw, na = na, dt = dz, amin = amin, amax = amax)
  Qrt_window[i,] <- ddj.t$sigma2.t
  D1t_window[i,] <- ddj.t$D1.t
  K4t_window[i,] <- ddj.t$K4.t
  D2Jumpt_window[i,] <- ddj.t$D2Jump.t
  lambdat_window[i,] <- ddj.t$lambda.t
  Jumpt_window[i,] <- ddj.t$Jump.t
}

zz = c(t(z_window))
R_tr = c(t(R_tr_window))
R_detr = c(t(R_detr_window))
Qrt = c(t(Qrt_window))
D1t = c(t(D1t_window))
K4t = c(t(K4t_window))
D2Jumpt = c(t(D2Jumpt_window))
lambdat = c(t(lambdat_window))
Jumpt = c(t(Jumpt_window))

ind1 <- which(z >= 700 & z < 740)

filename <- paste0("plots/snow_field_local_3_new.pdf")
pdf(filename,width = 6, height = 10)
par(mfrow=c(3,1),mar=c(5, 6, 2, 2),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(z[ind1],R[ind1], type = 'l',
     #ylim = c(-5,5),
     xlab=expression(italic(z)~"[mm]"),ylab=expression(italic(italic(R)^"'")), 
     #cex.lab = 1.3, cex.axis=1.3)
     cex.lab = 2, cex.axis=2)
plot(zz[ind1], K4t[ind1], type = 'l',
     ylim = c(0,1000),
     xlab=expression(italic(z)~"[mm]"),ylab=expression(italic(K)^"(4)"), 
     #cex.lab = 1.3, cex.axis=1.3)
     cex.lab = 2, cex.axis=2)
plot(zz[ind1],Qrt[ind1], type = 'l',
     ylim = c(0,2),
     #xlab=expression(italic(z)~"[mm]"),ylab=expression(italic(Q)), 
     xlab=expression(italic(z)~"[mm]"),ylab=expression(italic(sigma)[italic(xi)]^2), 
     #cex.lab = 1.3, cex.axis=1.3)
     cex.lab = 2, cex.axis=2)
dev.off()

###############################################################################################################################

##########################################################
#Snow Types Q criterion 
##########################################################
#rm(list = ls())
#graphics.off()

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
bw <- 0.3

# snow hardness data to analyze
data <- read.xls(xls ="data/Sample_profiles.xls")
z <- data$depth..mm.   # depth of the snow profile
#N <- length(z)  # length of the data
PP1 <- data$X2017.01.09.Sample07..kPa. # in table as PP
PP1o <- na.omit(PP1)
PP2 <- data$X2017.01.31.Sample06..kPa.
PP2o <- na.omit(PP2)
DH1 <- data$X2017.01.31.Sample05..kPa. #
DH1o <- na.omit(DH1)
DH2 <- data$X2017.01.31.Sample11..kPa. #
DH2o <- na.omit(DH2)
DH3 <- data$X2017.01.31.Sample10..kPa.
DH3o <- na.omit(DH3)
RG1 <- data$X2017.01.31.Sample01..kPa. #
RG1o <- na.omit(RG1)
RG2 <- data$X2017.01.31.Sample04..kPa. #
RG2o <- na.omit(RG2)
RG3 <- data$X2017.01.31.Sample03..kPa.
RG3o <- na.omit(RG3)
RGlr1 <- data$X2017.01.09.Sample04..kPa. #
RGlr1o <- na.omit(RGlr1)
RGlr2 <- data$X2017.01.09.Sample05..kPa. #
RGlr2o <- na.omit(RGlr2)
RGlr3 <- data$X2017.01.31.Sample07..kPa. #
RGlr3o <- na.omit(RGlr3)
RGlr4 <- data$X2017.01.09.Sample01..kPa.
RGlr4o <- na.omit(RGlr4)
RGlr5 <- data$X2017.01.09.Sample02..kPa.
RGlr5o <- na.omit(RGlr5)
RGlr6 <- data$X2017.01.09.Sample03..kPa.
RGlr6o <- na.omit(RGlr6)

######
# PP #
######

z1 <- z[1:length(PP1o)]
j <- which(z1>4)
zPP1 <- z1 <- z1[j]                  
PP1o <- PP1o[j]#/(1000*A)
PP1 <- PP1o

z2 <- z[1:length(PP2o)]
j <- which(z2>4)
zPP2 <- z2 <- z2[j]                  
PP2o <- PP2o[j]#/(1000*A)
PP2 <- PP2o

# calculate the trend of profile using Gaussian kernel
PP1_tr <- ksmooth(z1, PP1, "normal", bandwidth = 0.6)$y   
PP2_tr <- ksmooth(z2, PP2, "normal", bandwidth = 0.6)$y   

# detrend and normalize the snow profile
PP1 <- PP1 - PP1_tr
PP1 <- PP1/sd(PP1);
PP2 <- PP2 - PP2_tr
PP2 <- PP2/sd(PP2);

PP <- PP2
zz <- z2
steps <- c(1:25)#2^c(0:9)#1.5^c(0:15)#c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500)#c(1:10)
#l <- c(1:(length(PP1)-max(steps)))
N <- length(PP)
na <- 20   # number of points for computing the kernel
dz <- zz[2] - zz[1]    # depth step
amin <- min(PP)
amax <- max(PP)
avec <- seq(amin,amax,length.out=na)
dPP <- (amax-amin)/na

PPts <- cbind(zPP2,PP2)
ddj.t <- ddj_PP.t <- ddj(PPts,bandwidth = bw, na = na, dt = dz, amin = amin, amax = amax)


######
# DH #
######

z1 <- z[1:length(DH1o)]
j <- which(z1>4)
zDH1 <- z1 <- z1[j]                  
DH1o <- DH1o[j]#/(1000*A)
DH1 <- DH1o

z2 <- z[1:length(DH2o)]
j <- which(z2>4)
zDH2 <- z2 <- z2[j]                  
DH2o <- DH2o[j]#/(1000*A)
DH2 <- DH2o

z3 <- z[1:length(DH3o)]
j <- which(z3>4)
zDH3 <- z3 <- z3[j]                  
DH3o <- DH3o[j]#/(1000*A)
DH3 <- DH3o

# calculate the trend of profile using Gaussian kernel
DH1_tr <- ksmooth(z1, DH1, "normal", bandwidth = 0.6)$y   
DH2_tr <- ksmooth(z2, DH2, "normal", bandwidth = 0.6)$y
DH3_tr <- ksmooth(z3, DH3, "normal", bandwidth = 0.6)$y   

# detrend and normalize the snow profile
DH1 <- DH1 - DH1_tr
DH1 <- DH1/sd(DH1);
DH2 <- DH2 - DH2_tr
DH2 <- DH2/sd(DH2);
DH3 <- DH3 - DH3_tr
DH3 <- DH3/sd(DH3);

DH <- DH2
zz <- z2
steps <- c(1:25)#2^c(0:9)#1.5^c(0:15)#c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500)#c(1:10)
#l <- c(1:(length(DH1)-max(steps)))
N <- length(DH)
na <- 20   # number of points for computing the kernel
dz <- zz[2] - zz[1]    # depth step
amin <- min(DH)
amax <- max(DH)
avec <- seq(amin,amax,length.out=na)
dDH <- (amax-amin)/na

DHts <- cbind(zDH2,DH2)
ddj.t <- ddj_DH.t <- ddj(DHts,bandwidth = bw, na = na, dt = dz, amin = amin, amax = amax)


######
# RG #
######

z1 <- z[1:length(RG1o)]
j <- which(z1>4)
zRG1 <- z1 <- z1[j]                  
RG1o <- RG1o[j]#/(1000*A)
RG1 <- RG1o

z2 <- z[1:length(RG2o)]
j <- which(z2>4)
zRG2 <- z2 <- z2[j]                  
RG2o <- RG2o[j]#/(1000*A)
RG2 <- RG2o

z3 <- z[1:length(RG3o)]
j <- which(z3>4)
zRG3 <- z3 <- z3[j]                  
RG3o <- RG3o[j]#/(1000*A)
RG3 <- RG3o

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

RG <- RG2
zz <- z2
steps <- c(1:25)#2^c(0:9)#1.5^c(0:15)#c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500)#c(1:10)
#l <- c(1:(length(RG1)-max(steps)))
N <- length(RG)
na <- 20   # number of points for computing the kernel
dz <- zz[2] - zz[1]    # depth step
amin <- min(RG)
amax <- max(RG)
avec <- seq(amin,amax,length.out=na)
dRG <- (amax-amin)/na

RGts <- cbind(zRG2,RG2)
ddj.t <- ddj_RG.t <- ddj(RGts,bandwidth = bw, na = na, dt = dz, amin = amin, amax = amax)

########
# RGlr #
########

z1 <- z[1:length(RGlr1o)]
j <- which(z1>4)
zRGlr1 <- z1 <- z1[j]                  
RGlr1o <- RGlr1o[j]#/(1000*A)
RGlr1 <- RGlr1o

z2 <- z[1:length(RGlr2o)]
j <- which(z2>4)
zRGlr2 <- z2 <- z2[j]                  
RGlr2o <- RGlr2o[j]#/(1000*A)
RGlr2 <- RGlr2o

z3 <- z[1:length(RGlr3o)]
j <- which(z3>4)
zRGlr3 <- z3 <- z3[j]                  
RGlr3o <- RGlr3o[j]#/(1000*A)
RGlr3 <- RGlr3o

z4 <- z[1:length(RGlr4o)]
j <- which(z4>4)
zRGlr4 <- z4 <- z4[j]                  
RGlr4o <- RGlr4o[j]#/(1000*A)
RGlr4 <- RGlr4o

z5 <- z[1:length(RGlr5o)]
j <- which(z5>4)
zRGlr5 <- z5 <- z5[j]                  
RGlr5o <- RGlr5o[j]#/(1000*A)
RGlr5 <- RGlr5o

z6 <- z[1:length(RGlr6o)]
j <- which(z6>4)
zRGlr6 <- z6 <- z6[j]                  
RGlr6o <- RGlr6o[j]#/(1000*A)
RGlr6 <- RGlr6o

# calculate the trend of profile using Gaussian kernel
RGlr1_tr <- ksmooth(z1, RGlr1, "normal", bandwidth = 0.6)$y   
RGlr2_tr <- ksmooth(z2, RGlr2, "normal", bandwidth = 0.6)$y
RGlr3_tr <- ksmooth(z3, RGlr3, "normal", bandwidth = 0.6)$y   
RGlr4_tr <- ksmooth(z4, RGlr4, "normal", bandwidth = 0.6)$y   
RGlr5_tr <- ksmooth(z5, RGlr5, "normal", bandwidth = 0.6)$y   
RGlr6_tr <- ksmooth(z6, RGlr6, "normal", bandwidth = 0.6)$y   

# detrend and normalize the snow profile
RGlr1 <- RGlr1 - RGlr1_tr
RGlr1 <- RGlr1/sd(RGlr1);
RGlr2 <- RGlr2 - RGlr2_tr
RGlr2 <- RGlr2/sd(RGlr2);
RGlr3 <- RGlr3 - RGlr3_tr
RGlr3 <- RGlr3/sd(RGlr3);
RGlr4 <- RGlr4 - RGlr4_tr
RGlr4 <- RGlr4/sd(RGlr4);
RGlr5 <- RGlr5 - RGlr5_tr
RGlr5 <- RGlr5/sd(RGlr5);
RGlr6 <- RGlr6 - RGlr6_tr
RGlr6 <- RGlr6/sd(RGlr6);

RGlr <- RGlr5
zz <- z5
steps <- c(1:30)#2^c(0:9)#1.5^c(0:15)#c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500)#c(1:10)
#l <- c(1:(length(RGlr1)-max(steps)))
N <- length(RGlr)
na <- 20   # number of points for computing the kernel
dz <- zz[2] - zz[1]    # depth step
amin <- min(RGlr)
amax <- max(RGlr)
avec <- seq(amin,amax,length.out=na)
dRGlr <- (amax-amin)/na

RGlrts <- cbind(zRGlr5,RGlr5)
ddj.t <- ddj_RGlr.t <- ddj(RGlrts,bandwidth = bw, na = na, dt = dz, amin = amin, amax = amax)


#R <- c(PP2o, DH2o, RG2o, RGlr5o)
R <- c(PP2o,rep(NA,250))
R <- c(R, DH2o, rep(NA,311))
R <- c(R, RG2o, rep(NA,250))
R <- c(R, RGlr5o)
z <- c(0:(length(R)-1))*dz
K4.t <- c(ddj_PP.t$K4.t,rep(NA,250), ddj_DH.t$K4.t, rep(NA,311), ddj_RG.t$K4.t, rep(NA,250), ddj_RGlr.t$K4.t)
sigma2.t <- c(ddj_PP.t$sigma2.t,rep(NA,250), ddj_DH.t$sigma2.t, rep(NA,311), ddj_RG.t$sigma2.t, rep(NA,250), ddj_RGlr.t$sigma2.t)

filename <- paste0("plots/snow_types_local_3.pdf")
pdf(filename,width = 6, height = 10)
par(mfrow=c(3,1),mar=c(5, 6, 2, 2),mgp=c(3,1,0),oma=c(1,1,1,1))
plot(z,R, type = 'l',
     #log = 'xy', xaxt = "n", yaxt = "n", ylim = c(0.01,10),
     ylim = c(-5,250),
     xlab=expression(italic(z)~"[mm]"),ylab = expression(italic(R)~"[a.u.]"), 
     #cex.lab = 1.3, cex.axis=1.3)
     cex.lab = 2, cex.axis=2)
text(3, 210, expression(PP), col = 'red', cex = 2)
text(8.75, 210, expression(DH), col = 'red', cex = 2)
text(14.5, 210, expression(RG), col = 'red', cex = 2)
text(22, 210, expression(RGlr), col = 'red', cex = 2)
plot(z,K4.t, type = 'l',
     ylim = c(0,1000),
     xlab=expression(italic(z)~"[mm]"),ylab=expression(italic(K)^"(4)"), 
     #cex.lab = 1.3, cex.axis=1.3)
     cex.lab = 2, cex.axis=2)
plot(z, sigma2.t, type = 'l',
     ylim = c(0,2),
     #xlab=expression(italic(z)~"[mm]"),ylab=expression(italic(Q)), 
     xlab=expression(italic(z)~"[mm]"),ylab=expression(italic(sigma)[italic(xi)]^2), 
     #cex.lab = 1.3, cex.axis=1.3)
     cex.lab = 2, cex.axis=2)
dev.off()


