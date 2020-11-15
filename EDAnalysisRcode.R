############# DATA EXPLORATION AND PROCESSING ##########################
edvoldata = read.csv("EGDailyVolume.csv",header=T)
## Process Dates 
year = edvoldata$Year
month = edvoldata$Month
day = edvoldata$Day
datemat = cbind(as.character(day),as.character(month),as.character(year))
paste.dates = function(date){
    day = date[1]; month=date[2]; year = date[3]
    return(paste(day,month,year,sep="/"))
 }
dates = apply(datemat,1,paste.dates)
dates = as.Date(dates, format="%d/%m/%Y")
edvoldata = cbind(dates,edvoldata)
attach(edvoldata)

library(ggplot2)
ggplot(edvoldata, aes(dates, Volume)) + geom_line() + xlab("Time") + ylab("Daily ED Volume")

## ED Volume is count data: Transform
Volume.tr = sqrt(Volume+3/8)
hist(Volume,nclass=20,xlab="ED Volume", main="",col="brown")
hist(Volume.tr,nclass=20,xlab= "Transformed ED Volume", main="",col="blue")
ggplot(edvoldata, aes(dates, sqrt(Volume+3/8))) + geom_line() + xlab("Time") + ylab("Transformed Daily ED Volume")

################ TREND AND SEASONALITY ESTIMATION #########################
library(mgcv)

time.pts = c(1:length(Volume))
time.pts = c(time.pts - min(time.pts))/max(time.pts)
## Trend Estimation: Is there a trend?
## Local Polynomial Trend Estimation
loc.fit = loess(Volume.tr~time.pts)
vol.fit.loc = fitted(loc.fit)
## Splines Trend Estimation
gam.fit = gam(Volume.tr~s(time.pts))
summary(gam.fit)
vol.fit.gam = fitted(gam.fit)
## Is there a trend? 
ggplot(edvoldata, aes(dates, sqrt(Volume+3/8))) + geom_line() + xlab("Time") + ylab("Transformed Daily ED Volume")
lines(dates,vol.fit.loc,lwd=2,col="brown")
lines(dates,vol.fit.gam,lwd=2,col="red")


## Model Trend + Monthly Seasonality
## Using nonparametric trend and linear regression seasonality 
month = as.factor(format(dates,"%b"))
gam.fit.seastr.1 = gam(Volume.tr~s(time.pts)+month)
summary(gam.fit.seastr.1)
vol.fit.gam.seastr.1 = fitted(gam.fit.seastr.1)
ggplot(edvoldata, aes(dates, sqrt(Volume+3/8))) + geom_line() + xlab("Time") + ylab("Transformed Daily ED Volume")
lines(dates,vol.fit.gam.seastr.1,lwd=2,col="red")

## Add day-of-the-week seasonality
week = as.factor(weekdays(dates))
gam.fit.seastr.2 = gam(Volume.tr~s(time.pts)+month+week)
summary(gam.fit.seastr.2)
vol.fit.gam.seastr.2 = fitted(gam.fit.seastr.2)
## Compare the two fits: with & without day-of-the-week seasonality
ggplot(edvoldata, aes(dates, vol.fit.gam.seastr.2)) + geom_line() + xlab("Time") + ylab("Seasonality and Trend: Daily ED Volume")
lines(dates,vol.fit.gam.seastr.1,lwd=2,col="red")

## Does the addition of seasonality of day of the week adds predictive power?
lm.fit.seastr.1 = lm(Volume.tr~month)
lm.fit.seastr.2 = lm(Volume.tr~month+week)
anova(lm.fit.seastr.1,lm.fit.seastr.2)
vol.fit.lm.seastr.2 = fitted(lm.fit.seastr.2)
## Compare with & without trend
ggplot(edvoldata, aes(dates, vol.fit.gam.seastr.2)) + geom_line() + xlab("Time") + ylab("Seasonality and Trend: Daily ED Volume")
lines(dates,vol.fit.lm.seastr.2,lwd=2,col="blue")
lines(dates,vol.fit.gam,lwd=2,col="red")



################## STATIONARITY: RESIDUAL PROCESS ####################
## Residual Process: Trend Removal
resid.1 = Volume.tr-vol.fit.gam
## Residual Process: Stationarity Removal
resid.2 = Volume.tr-vol.fit.lm.seastr.2
## Residual Process: Trend & Stationarity Removal
resid.3 = Volume.tr-vol.fit.gam.seastr.2
y.min = min(c(resid.1,resid.2,resid.3))
y.max = max(c(resid.1,resid.2,resid.3))

ggplot(edvoldata, aes(dates, resid.1),ymin=y.min,ymax=y.max) + geom_line() + xlab("Time") + ylab("Residual Process")
lines(dates,resid.2,col="blue")
lines(dates,resid.3,col="brown")
#legend(2012,-3.5,legend=c("Trend","Season","Trend+Season"),lty = 1, col=c("black","blue","brown"))

acf(resid.1,lag.max=12*4,main="")
acf(resid.2,lag.max=12*4,main="",col="blue")
acf(resid.3,lag.max=12*4,main="",col="brown")
