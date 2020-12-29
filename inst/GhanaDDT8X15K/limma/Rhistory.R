#################################
######### PLOTS #################
#################################
par(mfrow=c(2,3), bg="lavender")#
thik = 2.5#
##############################
plot(S~time, data=LogSIRout, type="l", col="navy", lwd=thik, ylim=c(0,max(S)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~time, data=LogSIRout, type="l", col="darkred", lwd=thik, ylim=c(0,max(I)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(R~time, data=LogSIRout, type="l", col="darkgreen", lwd=thik, ylim=c(0,max(R)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(N~time, data=LogSIRout, type="l", col=5, lwd=thik, ylim=c(0,max(N)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(QOI~time, data=LogSIRout, type="l", col=6, lwd=thik, ylim=c(0,max(QOI)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~S, data=LogSIRout, type="l", col=1, lwd=thik)#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
#
#
###########################################
### Calculate Basic Reproductive Ratio ####
### and Equilibrium S-class ###############
###########################################
R.o <- unname(pars["beta"] / (pars["gamma"]+pars["mu"]+pars["nu"]))#
S.star <- 1/R.o#
R.o#
S.star#
#
#
#########################
### Equilibrium #########
### If S.star == final.S #
### We have Equilibrium#
#########################
Eqm <- LogSIRout[which(abs(LogSIRout[,"S"] - S.star) < 10e-6)[1],"time"]#
if (is.na(Eqm)) {cat("You have not yet reached Equilibrium!")#
} else {cat("Equilibrium is reached at t =", Eqm, "\n")}
##################
tInt <- seq(0, 75, by=1/4)#
pars <- c(r=0.5, beta=0.1, gamma=0.02, mu=0.2, nu=0.1, K=1000)#
initial <- c(S=5, I=5, R=5)#
#
LogSIR <- function(t, y, parms) { #
S = y[1] #
I = y[2]#
R = y[3]#
N = sum(y)#
QOI = y[2] / sum(y)       # set QOI here#
with(as.list(parms), {#
dS <- r*N*(1-N/K) - beta*S*I - mu*S#
dI <- beta*S*I  - gamma*I - mu*I - nu*I#
dR <- gamma*I - mu*R#
ODEs <- c(dS, dI, dR)#
list(ODEs, c(N, QOI))#
})#
}#
#
LogSIRout <- as.data.frame(lsoda(initial, times= tInt, parms= pars, func= LogSIR))#
colnames(LogSIRout)[5:6] <- c("N","QOI")   # rename 5th & 6th cols#
head(LogSIRout)#
#
#
#
#################################
######### PLOTS #################
#################################
par(mfrow=c(2,3), bg="lavender")#
thik = 2.5#
##############################
plot(S~time, data=LogSIRout, type="l", col="navy", lwd=thik, ylim=c(0,max(S)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~time, data=LogSIRout, type="l", col="darkred", lwd=thik, ylim=c(0,max(I)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(R~time, data=LogSIRout, type="l", col="darkgreen", lwd=thik, ylim=c(0,max(R)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(N~time, data=LogSIRout, type="l", col=5, lwd=thik, ylim=c(0,max(N)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(QOI~time, data=LogSIRout, type="l", col=6, lwd=thik, ylim=c(0,max(QOI)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~S, data=LogSIRout, type="l", col=1, lwd=thik)#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
#
#
###########################################
### Calculate Basic Reproductive Ratio ####
### and Equilibrium S-class ###############
###########################################
R.o <- unname(pars["beta"] / (pars["gamma"]+pars["mu"]+pars["nu"]))#
S.star <- 1/R.o#
R.o#
S.star#
#
#
#########################
### Equilibrium #########
### If S.star == final.S #
### We have Equilibrium#
#########################
Eqm <- LogSIRout[which(abs(LogSIRout[,"S"] - S.star) < 10e-6)[1],"time"]#
if (is.na(Eqm)) {cat("You have not yet reached Equilibrium!")#
} else {cat("Equilibrium is reached at t =", Eqm, "\n")}
LogSIRout
#
rm(list=ls())#
######################
require(deSolve)#
require(odesolve)
#
tInt <- seq(0, 75, by=1/4)#
pars <- c(r=0.5, beta=0.1, gamma=0.02, mu=0.2, nu=0.1, K=1000)#
initial <- c(S=5, I=5, R=5)#
#
LogSIR <- function(t, y, parms) { #
S = y[1] #
I = y[2]#
R = y[3]#
N = sum(y)#
QOI = y[2] / sum(y)       # set QOI here#
with(as.list(parms), {#
dS <- r*N*(1-N/K) - beta*S*I - mu*S#
dI <- beta*S*I  - gamma*I - mu*I - nu*I#
dR <- gamma*I - mu*R#
ODEs <- c(dS, dI, dR)#
list(ODEs, c(N, QOI))#
})#
}#
#
LogSIRout <- as.data.frame(lsoda(initial, times= tInt, parms= pars, func= LogSIR))#
colnames(LogSIRout)[5:6] <- c("N","QOI")   # rename 5th & 6th cols#
head(LogSIRout)#
#
#
#
#################################
######### PLOTS #################
#################################
par(mfrow=c(2,3), bg="lavender")#
thik = 2.5#
##############################
plot(S~time, data=LogSIRout, type="l", col="navy", lwd=thik, ylim=c(0,max(S)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~time, data=LogSIRout, type="l", col="darkred", lwd=thik, ylim=c(0,max(I)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(R~time, data=LogSIRout, type="l", col="darkgreen", lwd=thik, ylim=c(0,max(R)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(N~time, data=LogSIRout, type="l", col=5, lwd=thik, ylim=c(0,max(N)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(QOI~time, data=LogSIRout, type="l", col=6, lwd=thik, ylim=c(0,max(QOI)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~S, data=LogSIRout, type="l", col=1, lwd=thik)#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
#
#
###########################################
### Calculate Basic Reproductive Ratio ####
### and Equilibrium S-class ###############
###########################################
R.o <- unname(pars["beta"] / (pars["gamma"]+pars["mu"]+pars["nu"]))#
S.star <- 1/R.o#
R.o#
S.star#
#
#
#########################
### Equilibrium #########
### If S.star == final.S #
### We have Equilibrium#
#########################
Eqm <- LogSIRout[which(abs(LogSIRout[,"S"] - S.star) < 10e-6)[1],"time"]#
if (is.na(Eqm)) {cat("You have not yet reached Equilibrium!")#
} else {cat("Equilibrium is reached at t =", Eqm, "\n")}
################
tInt <- seq(0, 100, by=1/4)#
pars <- c(r=0.5, beta=0.1, gamma=0.02, mu=0.2, nu=0.1, K=1000)#
initial <- c(S=5, I=5, R=5)#
#
LogSIR <- function(t, y, parms) { #
S = y[1] #
I = y[2]#
R = y[3]#
N = sum(y)#
QOI = y[2] / sum(y)       # set QOI here#
with(as.list(parms), {#
dS <- r*N*(1-N/K) - beta*S*I - mu*S#
dI <- beta*S*I  - gamma*I - mu*I - nu*I#
dR <- gamma*I - mu*R#
ODEs <- c(dS, dI, dR)#
list(ODEs, c(N, QOI))#
})#
}#
#
LogSIRout <- as.data.frame(lsoda(initial, times= tInt, parms= pars, func= LogSIR))#
colnames(LogSIRout)[5:6] <- c("N","QOI")   # rename 5th & 6th cols#
head(LogSIRout)#
#
#
#
#################################
######### PLOTS #################
#################################
par(mfrow=c(2,3), bg="lavender")#
thik = 2.5#
##############################
plot(S~time, data=LogSIRout, type="l", col="navy", lwd=thik, ylim=c(0,max(S)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~time, data=LogSIRout, type="l", col="darkred", lwd=thik, ylim=c(0,max(I)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(R~time, data=LogSIRout, type="l", col="darkgreen", lwd=thik, ylim=c(0,max(R)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(N~time, data=LogSIRout, type="l", col=5, lwd=thik, ylim=c(0,max(N)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(QOI~time, data=LogSIRout, type="l", col=6, lwd=thik, ylim=c(0,max(QOI)))#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
plot(I~S, data=LogSIRout, type="l", col=1, lwd=thik)#
grid(NA, NULL,lty=2,lwd=0.5,col="gray90"); box()#
#
#
#
###########################################
### Calculate Basic Reproductive Ratio ####
### and Equilibrium S-class ###############
###########################################
R.o <- unname(pars["beta"] / (pars["gamma"]+pars["mu"]+pars["nu"]))#
S.star <- 1/R.o#
R.o#
S.star#
#
#
#########################
### Equilibrium #########
### If S.star == final.S #
### We have Equilibrium#
#########################
Eqm <- LogSIRout[which(abs(LogSIRout[,"S"] - S.star) < 10e-6)[1],"time"]#
if (is.na(Eqm)) {cat("You have not yet reached Equilibrium!")#
} else {cat("Equilibrium is reached at t =", Eqm, "\n")}
pine.data <- read.csv("Age vs DBH.csv", sep=",", header=T)
rm(list=ls())
pine.data <- read.csv("Age vs DBH.csv")
pine.data
fitPIFL <- lm(PIFL_dbh ~ age, data= pine.data)
summary(fitPIFL)     ### Summarized the model you just created; gives intercept, R-squared value, and P-value
anova.lm(fitPIFL)
fitPIAR <- lm(pine.data$PIAR_dbh ~ pine.data$age, data=pine.data)
summary(fitPIAR)
par(mfrow=c(1,2))#
plot(pine.data$PIFL_dbh ~ pine.data$age, data= pine.data); abline(fitPIFL, lwd=1.5)#
plot(pine.data$PIAR_dbh ~ pine.data$age, data= pine.data); abline(fitPIAR, lwd=1.5)
fitPIFL2 <- lm(PIFL_dbh ~ log(age), data= pine.data)#
fitPIAR2 <- lm(PIAR_dbh ~ log(age), data= pine.data)#
summary(fitPIFL2); summary(fitPIAR2)
PIFLm = coef(fitPIFL2)[2] # slope of fit2 PIFL#
PIFLb = coef(fitPIFL2)[1] # y-int of fit2 PIFL#
PIARm = coef(fitPIAR2)[2] # slope of fit2 PIAR#
PIARb = coef(fitPIAR2)[1] # y-int of fit2 PIAR#
PIFL.Res = residuals(fitPIFL2)#
PIAR.Res = residuals(fitPIAR2)#
PIFLm; PIFLb#
PIARm; PIARb
par(mfrow=c(1,2))
plot(pine.data$age[1:300], pine.data$PIFL_dbh[1:300], cex=0.75, xlab="Age", ylab="dbh PIFL")#
#scatter.smooth(pine.data$age, pine.data$PIFL.dbh, col=1)#
curve((PIFLm*log(x) + PIFLb), from = 0, to= 375, lty=4, add=T, col=2, lwd= 2)#
#
# Therefore PIFL dbh = 20.19726 * log(age) - 69.6496#
#
plot(pine.data$age, pine.data$PIAR_dbh, cex=0.75, xlab="Age", ylab="dbh PIAR")#
#scatter.smooth(pine.data$age, pine.data$PIFL.dbh, col=1)#
curve((PIARm*log(x) + PIARb), from = 0, to= 421, lty=4, add=T, col=2, lwd= 2)#
#
# Therefore PIAR dbh = 12.67591 * log(age) - 40.27601
x = c(0.1,0.15,0.22,0.34,0.5,0.75,1.1,1.7,2.6,3.8,5.8,8.9)#
y = c(-2.5,-2.2,-1.5,-1,-0.33,0,0.45,1.02,1.53,1.99,2.45,3.07)#
#
plot(x, y)#
#
fit <- lm(y ~ log(x))#
summary(fit)#
Res = summary(fit)$residuals#
m = coef(fit)[2] # slope#
b = coef(fit)[1] # y-int#
curve((m*log(x) + b), from = 0, to= 9, add=T, col=2)
nonlfit <- nls(y ~ (1-(b*exp(a-x))), start= list(a=0.1,b=1)) #can't seem to get good starting values#
summary(nonlfit)
nonlfit <- nls(y ~ (1-(b*exp(a-x))), start= list(a=0.1,b=1)) #can't seem to get good starting values
curve((1-(b*exp(a-x))), from=0, to= 9, lty= 4, col= 2, lwd= 2, add= T)
nonlfit <- nls(y ~ (1-(b*exp(a-x))), start= list(a=1,b=0.3))
m
b
fit.b <- coef(fit)[2]
plot(x, y, pch=19)
#
fit.a <- coef(fit)[1]#
fit.b <- coef(fit)[2]#
plot(x, y, pch=19)#
curve((1-(b*exp(a-x))), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
###############
x = seq(1, 26, 2); K= 0.4; m= 0.4; c= 15#
set.seed(125)#
y.det = K/(1+exp(m*(c-x)))#
#y <- jitter(y.det, amount=0.033)#
y <- rnorm(length(y.det), mean=y.det, sd=0.03)#
#logy = log(y)#
plot(x, y, main='Nonlinear Least Squares', xlab= 'time', ylab='Proportion Infected')
#################
fit <- nls(y ~ K/(1+exp(m*(c-x))), start= list(K=1, m=1, c=5))#
summary(fit)#
#
K1 <- coef(fit)["K"]#
m1 <- coef(fit)["m"]#
c1 <- coef(fit)["c"]
#
plot(x, y, main='Nonlinear Least Squares', pch=19, xlab= 'time', ylab='Proportion Infected')#
curve(K1/(1+exp(m1*(c1-x))), from=0, to=max(x), lty= 1, col=4, lwd=2, add=TRUE)#
curve(CI95[1,1]/(1+exp(CI95[2,2]*(CI95[3,2]-x))), from=0, to=max(x), lty= 4, col=2, lwd=2, add=TRUE)#
curve(CI95[1,2]/(1+exp(CI95[2,1]*(CI95[3,1]-x))), from=0, to=max(x), lty= 4, col=2, lwd=2, add=TRUE)
#################
fit <- nls(y ~ K/(1+exp(m*(c-x))), start= list(K=1, m=1, c=5))#
summary(fit)#
#
K1 <- coef(fit)["K"]#
m1 <- coef(fit)["m"]#
c1 <- coef(fit)["c"]#
#
CI95 <- confint(fit)#
#
plot(x, y, main='Nonlinear Least Squares', pch=19, xlab= 'time', ylab='Proportion Infected')#
curve(K1/(1+exp(m1*(c1-x))), from=0, to=max(x), lty= 1, col=4, lwd=2, add=TRUE)#
curve(CI95[1,1]/(1+exp(CI95[2,2]*(CI95[3,2]-x))), from=0, to=max(x), lty= 4, col=2, lwd=2, add=TRUE)#
curve(CI95[1,2]/(1+exp(CI95[2,1]*(CI95[3,1]-x))), from=0, to=max(x), lty= 4, col=2, lwd=2, add=TRUE)
fit <- nls(y ~ K/(1+exp(m*(c-x))), start= list(K=1, m=1, c=5))
x = c(0.1,0.15,0.22,0.34,0.5,0.75,1.1,1.7,2.6,3.8,5.8,8.9)#
y = c(-2.5,-2.2,-1.5,-1,-0.33,0,0.45,1.02,1.53,1.99,2.45,3.07)
y
x
fit <- lm(y ~ log(x))
summary(fit)
plot(y~x, pch=19)
fit.b
curve((1-(b*exp(a-x))), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
plot(y~x, pch=19)#
curve((1-(b*exp(a-x))), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
#
plot(x, y)#
#
fit <- lm(y ~ log(x))#
summary(fit)#
Res <- summary(fit)$residuals#
m <- coef(fit)[2] # slope#
b <- coef(fit)[1] # y-int#
curve((m*log(x) + b), from = 0, to= 9, add=TRUE, col=2)
curve(1 - (b*exp(a-x)), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
curve(1 - (b*exp(m-x)), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
plot(y~x, pch=19)#
curve(1 - (b*exp(m-x)), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
nonlfit <- nls(y ~ (1-(b*exp(a-x))), start= list(a=1, b=0.3))
nonlfit <- nls(y ~ 1-(b*exp(a-x)), start= list(a=1, b=0.3))
nonlfit <- nls(y ~ 1-(b*exp(a-x)), start= list(a=0.5, b=1))
nonlfit <- nls(y ~ 1-(b*exp(a-x)), start= list(a=0.5, b=2))
?try
nonlfit <- nls(y ~ try(1-(b*exp(a-x))), start= list(a=0.5, b=0.2))
nonlfit <- nls(y ~ 1-(b*exp(a-x)), start= list(a=0.5, b=0.2))
mydata <- data.frame(x,y)
mydata
plot(mydata)
nonlfit <- nls(y ~ 1-(b*exp(a-x)), start= list(a=0.5, b=0.2), data=mydata)
nonlfit <- nls(y ~ 1-(b*exp(a-x)), start= list(a=1, b=0.2), data=mydata)
nonlfit <- nls(y ~ 1-(exp(a-x)), start= list(a=1, b=0.2), data=mydata)
nonlfit <- nls(y ~ 1-(exp(a-x)), start= list(a=1), data=mydata)
summary(nonlfit)
fit.a <- coef(fit)[1]
fit.a
coef(nonlfit)
nonlfit.a <- coef(nonlfit)
a <- coef(nonlfit)#
plot(y~x, pch=19)#
curve(1 - (m*exp(b-x)), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)#
curve(1 - (exp(a-x)), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
plot(y~x, pch=19)#
curve(1 - (m*exp(b-x)), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)#
curve(1 - (exp(a-x)), from=0, to= 9, lty= 4, col= 2, lwd= 2, add=TRUE)
a <- coef(nonlfit)
a
plot(pine.data$age[1:300], PIFL.Res, cex=0.75, xlab="log(Age)", ylab="Residuals"); abline(h=0)#
plot(pine.data$age, PIAR.Res, cex=0.75, xlab="log(Age)", ylab="Residuals"); abline(h=0)
plot(pine.data$age, pine.data$PIAR_dbh, cex=0.75, xlab="Age", ylab="dbh PIAR")#
#scatter.smooth(pine.data$age, pine.data$PIFL.dbh, col=1)#
curve((PIARm*log(x) + PIARb), from = 0, to= 421, lty=4, add=T, col=2, lwd= 2)
plot(pine.data$age[1:300], pine.data$PIFL_dbh[1:300], cex=0.75, xlab="Age", ylab="dbh PIFL")#
#scatter.smooth(pine.data$age, pine.data$PIFL.dbh, col=1)#
curve((PIFLm*log(x) + PIFLb), from = 0, to= 375, lty=4, add=T, col=2, lwd= 2)
fitPIFL2 <- lm(PIFL_dbh ~ log(age), data= pine.data)#
fitPIAR2 <- lm(PIAR_dbh ~ log(age), data= pine.data)#
summary(fitPIFL2); summary(fitPIAR2)#
#
PIFLm = coef(fitPIFL2)[2] # slope of fit2 PIFL#
PIFLb = coef(fitPIFL2)[1] # y-int of fit2 PIFL#
PIARm = coef(fitPIAR2)[2] # slope of fit2 PIAR#
PIARb = coef(fitPIAR2)[1] # y-int of fit2 PIAR#
PIFL.Res = residuals(fitPIFL2)#
PIAR.Res = residuals(fitPIAR2)#
PIFLm; PIFLb#
PIARm; PIARb
plot(pine.data$PIFL_dbh ~ pine.data$age, data= pine.data); abline(fitPIFL, lwd=1.5)#
plot(pine.data$PIAR_dbh ~ pine.data$age, data= pine.data); abline(fitPIAR, lwd=1.5)#
#### NOT TOO GOOD A FIT - TRANSFORM
plot(pine.data$age, pine.data$PIFL_dbh)#
points(pine.data$age, pine.data$PIAR_dbh)#
curve((mean(c(20.19726,12.67591)))*log(x) + mean(c(-69.6496,-40.27601)), from = 0, to= 425, add=TRUE, col=2)
#########
rm(list=ls())#
require(rgl)#
require(emdbook)
x <- sort(rnorm(1000))#
y <- rnorm(1000)#
z <- rnorm(1000) + atan2(x,y)#
plot3d(x, y, z, col=rainbow(1000))
x <- rnorm(1000)
x <- rnorm(1000)#
y <- rnorm(1000)
z <- rnorm(1000)
x <- sort(rnorm(1000))
plot3d(x, y, z, col=rainbow(1000))
x <- rnorm(1000)#
y <- rnorm(1000)#
z <- rnorm(1000)#
plot3d(x, y, z, col=rainbow(1000))
#
x <- rnorm(1000)#
y <- rnorm(1000)#
z <- rnorm(1000)#
plot3d(sort(x), y, z, col=rainbow(1000))#
plot3d(x, y, z, col=rainbow(1000), type="s", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="s", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="shade", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="p", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="wire", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="dots", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="wireframe", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="l", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="h", size=0.7)
plot3d(sort(x), y, z, col=rainbow(1000), type="h", size=0.1)
?plot3d
plot3d(sort(x), y, z, col=rainbow(1000),add=T)
plot3d(sort(x), y, z, col=rainbow(1000))
plot3d(sort(x), y, z, col=rainbow(1000), type="h", size=0.1, add=T)
######
x <- rnorm(10000)#
y <- rnorm(10000)#
z <- rnorm(10000)#
plot3d(sort(x), y, z, col=rainbow(1000))#
plot3d(sort(x), y, z, col=rainbow(10000), type="s", size=0.1, add=T)
#
x <- rnorm(5000)#
y <- rnorm(5000)#
z <- rnorm(5000)#
plot3d(sort(x), y, z, col=rainbow(5000))
x <- rnorm(5000)#
y <- rnorm(5000)#
z <- rnorm(5000)#
plot3d(sort(x), y, z, col=rainbow(5000))
require(plot3d)
require(R.basic)
require(R.classic)
require(R.oo)
install.packages("R.classes", contriburl="http://www.maths.lth.se/help/R/R.classes/")
library(R.oo)
require(scatterplot3d)
scatterplot3d(sort(x), y, z, col=rainbow(5000))
library(Rcmdr)
attach(mtcars)
scatterplot3d(sort(x), y, z)
plot3d(sort(x), y, z, col=rainbow(5000))
plot3d(sort(x), y, z, col=rainbow(5000), type="s", size=0.5)
scatter3d(wt, disp, mpg)
adjPval
setwd("~/Documents/Dropbox/Black & Donnelly/LIMMA")
setwd("~/Documents/Dropbox/CSU/Black & Donnelly/LIMMA")
####################
require(limma)#
require(maanova)#
require(qvalue)#
require(calibrate)#
require(lattice)
source(paste(getwd(),"array.frame().R",sep="/")  ### <- set path here
)
ls()
source(paste(getwd(),"search.replace().R",sep="/"))   # call in search.replace()
