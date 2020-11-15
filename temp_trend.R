t_d = read.table("atl_temp.txt",header=T)
temp = as.vector(temp(data[,-c(1,14)]))
temp = ts(temp,start=1879,frequency=12)
head(temp)
ts.plot(temp,ylab="Temperature")

############### TREND ESTIMATION ###################
## Is there a trend in the average temperature?
time.1 = c(1:length(temp))
tail(time.1)
#min(time.1)/max(time.1) = 0.0006038647
time.pts = c(time.pts - min(time.pts))/max(time.pts)
tail(time.pts, 25)
## Fit a moving average 
mav.fit = ksmooth(time.pts, temp, kernel = "box")
temp.fit.mav = ts(mav.fit$y,start=1879,frequency=12)
ts.plot(temp,ylab="Temperature")
lines(temp.fit.mav,lwd=2,col="purple")
abline(temp.fit.mav[1],0,lwd=2,col="blue")

## Fit a parametric quadraric polynomial
x1 = time.pts
x2 = time.pts^2
lm.fit = lm(temp~x1+x2)
summary(lm.fit)
temp.fit.lm = ts(fitted(lm.fit),start=1879,frequency=12)
ts.plot(temp,ylab="Temperature")
lines(temp.fit.lm,lwd=2,col="green")
abline(temp.fit.mav[1],0,lwd=2,col="blue")

## Fit a trend using non-parametric regression
## Local Polynomial Trend Estimation
loc.fit = loess(temp~time.pts)
temp.fit.loc = ts(fitted(loc.fit),start=1879,frequency=12)
## Splines Trend Estimation
library(mgcv)
gam.fit = gam(temp~s(time.pts))
temp.fit.gam = ts(fitted(gam.fit),start=1879,frequency=12)
## Is there a trend? 
ts.plot(temp,ylab="Temperature")
lines(temp.fit.loc,lwd=2,col="brown")
lines(temp.fit.gam,lwd=2,col="red")
abline(temp.fit.loc[1],0,lwd=2,col="blue")

## Compare all estimated trends
all.val = c(temp.fit.mav,temp.fit.lm,temp.fit.gam,temp.fit.loc)
ylim= c(min(all.val),max(all.val))
ts.plot(temp.fit.lm,lwd=2,col="green",ylim=ylim,ylab="Temperature")
lines(temp.fit.mav,lwd=2,col="purple")
lines(temp.fit.gam,lwd=2,col="red")
lines(temp.fit.loc,lwd=2,col="brown")
legend(x=1900,y=64,legend=c("MAV","LM","GAM","LOESS"),lty = 1, col=c("purple","green","red","brown"))

################ SEASONALITY ESTIMATION #########################

library(dynlm)

## Estimate seasonality using ANOVA approach
month = season(temp)
## Drop January (model with intercept)
model1 = lm(temp~month)
summary(model1)
## All seasonal mean effects (model without intercept)
model2 = lm(temp~month-1)
summary(model2)

## Estimate seasonality using cos-sin model
har=dynlm(temp,1)
model3=lm(temp~har)
summary(model3)
har2=harmonic(temp,2)
model4=lm(temp~har2)
summary(model4)

## Compare Seasonality Estimates
## Seasonal Means Model
st1 = coef(model2)
## Cos-Sin Model
st2 = fitted(model4)[1:12]
plot(1:12,st1,lwd=2,type="l",xlab="Month",ylab="Seasonality")
lines(1:12,st2,lwd=2, col="brown")

################ TREND AND SEASONALITY ESTIMATION #########################
## Using linear regression

## Fit a parametric model for both trend and seasonality
x1 = time.pts
x2 = time.pts^2
har2=harmonic(temp,2)
lm.fit = lm(temp~x1+x2+har2)
summary(lm.fit)
dif.fit.lm = ts((temp-fitted(lm.fit)),start=1879,frequency=12)
ts.plot(dif.fit.lm,ylab="Residual Process")

## Fit a non-parametric model for trend and linear model for seasonality
gam.fit = gam(temp~s(time.pts)+har2)
dif.fit.gam = ts((temp-fitted(gam.fit)),start=1879,frequency=12)
ts.plot(dif.fit.gam,ylab="Residual Process")

## Compare approaches 
ts.plot(dif.fit.lm,ylab="Residual Process",col="brown")
lines(dif.fit.gam,col="blue")

acf(temp,lag.max=12*4,main="")
acf(dif.fit.lm,lag.max=12*4,main="")
acf(dif.fit.gam,lag.max=12*4,main="")

