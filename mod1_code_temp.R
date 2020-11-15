## A note on TSA -- this library is not available for more recent R versions
## We will use the following R founction from TSA -- creates the harmonic functions for the cos-sine seasonality model
harmonic <-
  function (x, m = 1) 
  {
    if (!is.ts(x) || (2 * m) > frequency(x)) 
      stop("x need to be a time series with 2m <= frequency(x)")
    y = outer(2 * pi * time(x), 1:m)
    cosy = apply(y, 2, cos)
    siny = apply(y, 2, sin)
    mult = 2 * (1:m)
    colnames(cosy) = paste(paste("cos(", mult, sep = ""), "*pi*t)", 
                           sep = "")
    colnames(siny) = paste(paste("sin(", mult, sep = ""), "*pi*t)", 
                           sep = "")
    out = cbind(cosy, siny)
    colnames(out) = c(colnames(cosy), colnames(siny))
    if ((2 * m) == frequency(x)) 
      out = out[, -(2 * m)]
    invisible(out)
  }


###################################################################################################################3
data = read.table("AvTempAtlanta.txt",header=T)
temp = as.vector(t(data[,-c(1,14)]))
temp = ts(temp,start=1879,frequency=12)
ts.plot(temp,ylab="Temperature")

############### TREND ESTIMATION ###################
## Is there a trend in the average temperature?

## X-axis points converted to 0-1 scale, common in nonparametric regression
time.pts = c(1:length(temp))
time.pts = c(time.pts - min(time.pts))/max(time.pts)

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
## Drop January (model with intercept)
model1 = dynlm(temp~season(temp))
summary(model1)
## All seasonal mean effects (model without intercept)
model2 = dynlm(temp~season(temp)-1)
summary(model2)

## Estimate seasonality using cos-sin model
model3=dynlm(temp~harmon(temp))
summary(model3)
model4=dynlm(temp~harmon(temp,2))
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
## Linear trend
lm.fit.lin = dynlm(temp~trend(temp)+harmon(temp,2))
## Quadratic trend
x1 = time.pts
x2 = time.pts^2
lm.fit = dynlm(temp~x1+x2+harmon(temp,2))
summary(lm.fit)
dif.fit.lm = ts((temp-fitted(lm.fit)),start=1879,frequency=12)
ts.plot(dif.fit.lm,ylab="Residual Process")

## Fit a non-parametric model for trend and linear model for seasonality
har2 = harmonic(temp,2)
gam.fit = gam(temp~s(time.pts)+har2)
dif.fit.gam = ts((temp-fitted(gam.fit)),start=1879,frequency=12)
ts.plot(dif.fit.gam,ylab="Residual Process")

## Compare approaches 
ts.plot(dif.fit.lm,ylab="Residual Process",col="brown")
lines(dif.fit.gam,col="blue")

acf(temp,lag.max=12*4,main="")
acf(dif.fit.lm,lag.max=12*4,main="")
acf(dif.fit.gam,lag.max=12*4,main="")

