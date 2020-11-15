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
Volume.tr = sqrt(Volume+3/8)

## Model Trend + Monthly Seasonality
library(mgcv)
time.pts = c(1:length(Volume))
time.pts = c(time.pts - min(time.pts))/max(time.pts)
month = as.factor(format(dates,"%b"))
week = as.factor(weekdays(dates))
gam.fit.seastr = gam(Volume.tr~s(time.pts)+month+week)
vol.fit.gam.seastr = fitted(gam.fit.seastr)
resid.process = Volume.tr-vol.fit.gam.seastr

################################################################
acf(resid.process,lag.max=12*4,main="ACF: Residual Plot")
pacf(resid.process,lag.max=12*4,main="PACF: Residual Plot")
###############################################################
library(TSA)
###### Fit an AR(p) process for for p<= order.max ####
mod = ar(resid.process,order.max=20)
# What is the selected order?
print(mod$order)
# find the list of arguments provided by AR fit
summary(mod)
# plot aic values
# On non-log scale it is difficult to detect the minimum.
plot(c(0:20),mod$aic+.0001, type="b",log="y",xlab="order",ylab="log-AIC values")

# Are the roots of fitted AR within the unit circle?
# extract roots from the model output
roots = polyroot(c(1,(-mod$ar)))
# adjust the x and y -axis limits to include the full circle
plot(roots,xlim=c(-2,1.5),ylim=c(-1.5,1.5))
# draw a unit root circle
lines(complex(arg = seq(0,2*pi,len=300)))
 
# residuals analysis
resids = mod$resid[(mod$order+1): length(mod$resid)]

par(mfrow=c(2,2))
plot(resids,xlab="",ylab="Model Residuals")
acf(resids,main='ACF of the Model Residuals')
pacf(resids,main='PACF of the Model Residuals')
qqnorm(resids)

###### Fit an ARMA(p,q) for some values of p and q #####
modarma = arima(resid.process, order = c(6,0,1),method = "ML")

par (mfrow=c(2,2))
plot(resid(modarma), ylab='Standardized Residuals')
abline(h=0)
acf(as.vector(resid(modarma)),main= 'ACF of the Model Residuals')
pacf(as.vector(resid(modarma)),main='PACF of the Model Residuals')
qqnorm(resid(modarma))
qqline(resid(modarma))

## Order selection -- EACF

eacf(resid.process,ar.max = 6, ma.max = 6)

## Order selection -- AIC 
n = length(resid.process)
norder = 6
p = c(1:norder)-1; q = c(1:norder)-1
aic = matrix(0,norder,norder)
for(i in 1:norder){
   for(j in 1:norder){
    modij = arima(resid.process,order = c(p[i],0,q[j]), method='ML')
    aic[i,j] = modij$aic-2*(p[i]+q[j]+1)+2*(p[i]+q[j]+1)*n/(n-p[i]-q[j]-2)
   }  
 }

aicv = as.vector(aic)  
plot(aicv,ylab="AIC values")
indexp = rep(c(1:norder),norder)
indexq = rep(c(1:norder),each=norder)
indexaic = which(aicv == min(aicv))
porder = indexp[indexaic]-1
qorder = indexq[indexaic]-1
final_model = arima(resid.process,order = c(porder,0,qorder), method='ML')

par (mfrow=c(2,2))
plot(resid(final_model), ylab='Standardized Residuals')
abline(h=0)
acf(as.vector(resid(final_model)),main= 'ACF of the Model Residuals')
pacf(as.vector(resid(final_model)),main='PACF of the Model Residuals')
qqnorm(resid(final_model))
qqline(resid(final_model))

#### Test for Independence for final model
Box.test(final_model$resid, lag = (porder+qorder+1), type = "Box-Pierce", fitdf = (porder+qorder))
Box.test(final_model$resid, lag = (porder+qorder+1), type = "Ljung-Box", fitdf = (porder+qorder))

#### Test for Independence for smaller model
Box.test(modarma$resid, lag = (porder+qorder+1), type = "Box-Pierce", fitdf = (porder+qorder))
Box.test(modarma$resid, lag = (porder+qorder+1), type = "Ljung-Box", fitdf = (porder+qorder))
