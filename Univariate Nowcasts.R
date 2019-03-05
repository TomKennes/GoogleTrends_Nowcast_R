###################################################################################
# CONSUMER CONFIDENCE USING GOOGLE TRENDS
# Author: Tom Kennes              
# Date: 03/12/2017           
###################################################################################

rm(list = ls())
cat("\014")
dates_creator <- function(x, y, s = T){
  n = (y[1] - x[1])*12 - (x[2] - y[2]) + 1
  dats <- t(as.matrix(c(as.character(x[1]), paste("0", as.character(x[2]), sep = ""))))
  i = 2
  while(i <= n){
    if(i == 2){
      last = dats
    } else{
      last = dats[nrow(dats),]
    }
    if(last[2] == "12"){
      new = c(as.character(as.numeric(last[1]) + 1), "01")
    } else {
      if(as.numeric(last[2]) >= 9){
        new = c(last[1], as.character(as.numeric(last[2]) + 1))
      } else {
        new = c(last[1], paste("0", as.character(as.numeric(last[2]) + 1), sep = ""))
      }
    }
    dats = rbind(dats, new)    
    i = i + 1
  }
  
  
  if(s){
    dates = c()
    for(i in 1:nrow(dats)){
      dates = c(dates, paste(paste(dats[i,1], dats[i,2], sep = "-"), "-01", sep = ""))
    }
    dates = as.Date(dates)
    return(dates)
  } else {
    colnames(dats) = c("years", "months")
    return(dats)
  }
}
## READ AND CLEAN DATA

#Clear workspace memory
set.seed = 1234

#Set working directory and load data
setwd("#######################################")
rawX <- read.csv2("GoogleQueries01.csv", check.names = FALSE)
rawy <- read.csv("CBS_cc.csv", header = F, sep = ",")

#Set sizes
n <- nrow(rawX)
p <- ncol(rawX)
#First two rows of X are irrelevant
rawX <- rawX[,3:p]
rawy <- as.matrix(rawy)
p <- p-2

library("tseries")
#Test for unit roots and if so, take the difference
if(adf.test(rawy)$p.value > 0.01)
{
  y <- diff(rawy)
}
X <- diff(as.matrix(rawX))
for(i in 1:p)
{
  if(adf.test(rawX[,i])$p.value <= 0.01)
  {
    #process is stationary, replace in X
    X[,i] <- rawX[2:n,i]
    print("stationary")
  }
}
n <- nrow(X)

#Select order of autoregression and include in the set of explanatory variables
ar(y)$order
lag.y <- y[1:n-1]
y <- y[1:n]
X <- X[2:n,]
X <- cbind(lag.y,X)
n <- nrow(X)
p <- ncol(X)

#Split sample into estimation and validation
Xest <- X[1:94,]
yest <- y[1:94]
Xval <- X[95:n,]
yval <- y[95:n]
yest <- as.matrix(yest)

#Standardize the variables
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
Xmeanest <- colMeans(Xest)
Xsdest <- colSd(Xest)
ymeanest <- mean(yest)
ysdest <- sd(yest)
Xest <- scale(Xest, TRUE, TRUE)
yest <- scale(yest, TRUE, TRUE)
Xval <- scale(Xval,center = Xmeanest,scale = Xsdest)
yval <- scale(yval,center = ymeanest,scale = ysdest)
cc = scale(y, center = ymeanest, scale = ysdest)
cc_test = cc[96:length(cc)]
cc_training = cc[1:95]

# Univariate Forecasting
x= c(2008,1)
y = c(2017, 10)
dates = dates_creator(x, y)
data = data.frame(cc = as.vector(cc), dates = dates)
require(ggplot2)
ggplot( data = data, aes(x = dates, y = cc)) + 
  geom_line(size = 2, color = "#0072B2") + 
  ggtitle("Consumer Confidence") + 
  xlab("Time") + 
  ylab("Value")
data = data[96:118,]

require(ggplot2)
ggplot( data = data, aes(x = dates, y = cc)) + 
  geom_line(size = 2, color = "#0072B2") + 
  ggtitle("Consumer Confidence") + 
  xlab("Time") + 
  ylab("Value")

# Time series modelling
require(stats)

lenq = 5
lenp = 5
lend = 2

spec <- rep(0, (lenq + 1)*(lenp + 1)*(lend + 1))
NN = length(spec)
loglike <- spec
aic <- spec
bic <- spec
hq <- spec
count <- 1
for(p in 0:lenp){
  for(q in 0:lenq){
    for(d in 0:lend){
      #THE LAST ENTRY OF arima(INPUT[...,x])
      a <- arima(cc_training, order = c(p,d,q))
      loglike[count] <- a$loglik
      aic[count] <- a$aic
      bic[count] <- -2*a$loglik + 2*log((2*a$loglik + a$aic)/2)
      hq[count] <-  -2*a$loglik + 2*log(log((2*a$loglik + a$aic)/2))
      spec[count] <- paste(paste(as.character(p), as.character(d)), as.character(q))
      count <- count + 1
      print(count/NN)
    }
  }
}

avg_ic <- c()
for(i in 1:NN){
  avg_ic <- c(avg_ic, (aic[i] + bic[i] + hq[i])/3)
}
print(avg_ic)

order_ = sort(bic, decreasing = FALSE, index.return = T)$ix

best_models <- cbind(spec[order_],
                     avg_ic[order_],
                     aic [order_],
                     bic [order_],
                     hq  [order_],
                     loglike [order_])

best_models <- best_models[best_models[,2] != -Inf,]
best_models <- best_models[best_models[,2] != -Inf,]
colnames(best_models) <- c("ARIMA", "avg IC", "AIC", "BIC", "HQ", "Loglik")
print(best_models[1:10,])
continue = best_models[,1]

msfe_ <- c()

forecasts <- cc_test
coln <- c("Actual")
require(forecast)
pdf("All forecasts.pdf")
for(i in 1:length(continue)){
  spec = as.numeric(strsplit(continue[i], " ")[[1]])
  spec_ = paste(spec, collapse = ", ")
  coln <- c(coln, spec_)
  use = arima(cc_training, order = spec)
  fcast = forecast(use, h = length(cc_test))
  ff_ <- cbind(cc_training, fcast$fitted, rep(NA, length(cc_training)), rep(NA, length(cc_training)), rep(NA, length(cc_training)), rep(NA, length(cc_training)))
  ff <- cbind(cc_test, fcast$mean, fcast$lower, fcast$upper)
  cc__ = colnames(ff)
  colnames(ff_) <- cc__
  datatjes <- rbind(ff_, ff)
  datatjes <- data.frame(cbind(dates, datatjes))
  datatjes[,1] = dates
  colnames(datatjes)[2] = "cc"
  forecasts = cbind(forecasts, fcast$mean)
  fmse = mean((cc_training - fcast$fitted)^2)
  msfe = mean((cc_test - as.vector(fcast$mean))^2)
  msfe_ <- c(msfe_, msfe)
  print(ggplot(data = datatjes, aes(x = dates, y = cc))+
    geom_line(col = "red") +
    geom_line(aes(y = fcast.mean), col = "blue") +
    geom_ribbon(aes(ymin = fcast.lower.95., ymax = fcast.upper.95.), alpha = 0.25) +
    ggtitle(paste("Consumer Confidence. Fitted and forecasted with \n Arima(", spec_, ")", sep = "")) +
    xlab(paste("Dates.\t fitted MSE: ", as.character(round(fmse,3))," MSFE: ", as.character(round(msfe, 3)), sep = "")) + 
    ylab("values"))
}
dev.off()



require(forecast)
take <- sort(msfe_, decreasing = F, index.return = T)$ix[1:10]
continue = continue[take]
continue
pdf("Time Series Forecast Benchmark.pdf")
print(ggplot( data = data, aes( dates, cc )) + 
  geom_line(size = 2, color = "#0072B2") + 
  ggtitle("Consumer Confidence") + 
  xlab("Time") + 
  ylab("Value"))
for(i in 1:length(continue)){
  spec = as.numeric(strsplit(continue[i], " ")[[1]])
  spec_ = paste(spec, collapse = ", ")
  use = arima(cc_training, order = spec)
  fcast = forecast(use, h = length(cc_test))
  ff_ <- cbind(cc_training, fcast$fitted, rep(NA, length(cc_training)), rep(NA, length(cc_training)), rep(NA, length(cc_training)), rep(NA, length(cc_training)))
  ff <- cbind(cc_test, fcast$mean, fcast$lower, fcast$upper)
  cc__ = colnames(ff)
  colnames(ff_) <- cc__
  datatjes <- rbind(ff_, ff)
  datatjes <- data.frame(cbind(dates, datatjes))
  datatjes[,1] = dates
  colnames(datatjes)[2] = "cc"
  fmse = mean((cc_training - fcast$fitted)^2)
  msfe = mean((cc_test - as.vector(fcast$mean))^2)
  msfe_ <- c(msfe_, msfe)
  print(ggplot(data = datatjes, aes(x = dates, y = cc))+
          geom_line(col = "red") +
          geom_line(aes(y = fcast.mean), col = "blue") +
          geom_ribbon(aes(ymin = fcast.lower.95., ymax = fcast.upper.95.), alpha = 0.25) +
          ggtitle(paste("Consumer Confidence. Fitted and forecasted with \n Arima(", spec_, ")", sep = "")) +
          xlab(paste("Dates.\t fitted MSE: ", as.character(round(fmse,3))," E[MSFE]: ", as.character(round(msfe, 3)), 
                     "\n mean Interval: [",as.character(round(mean(fcast$lower[,2]),3)),", ", as.character(round(mean(fcast$upper[,2]),3)), "]",  sep = "")) + 
          ylab("values"))
}
dev.off()

colnames(forecasts) = coln
write.csv(forecasts, "Forecasts Univariate.csv")
pick = 6
take.first.five = T
pdf("Forecasts 6step.pdf")
sets = c(1,0,0)
for(i in 1:3){
  sets = rbind(sets, as.numeric(strsplit(colnames(forecasts)[i + 1], ", ")[[1]]))
}
for(set in 1:nrow(sets)){
  tmp <- c(0,0,0,0,0)
  for(i in 0:(length(cc_test) - pick)){
    if(i > 0){
      train = c(cc_training, cc_test[1:i])
    } else {
      train = cc_training
    }
    use = arima(train, order = sets[set,])
    fore = forecast(use, h = pick)
    if(i == 0 && take.first.five){
      tmp = rbind(tmp, cbind(fore$mean, fore$lower, fore$upper))
    } else {
      tmp = rbind(tmp, c(fore$mean[pick], fore$lower[pick,], fore$upper[pick,]))
    }
  }
  tmp <- tmp[2:nrow(tmp), ]
  ff_ <- cbind(cc_training, cc_training, rep(NA, length(cc_training)), rep(NA, length(cc_training)), rep(NA, length(cc_training)), rep(NA, length(cc_training)))
  ff <- cbind(cc_test, tmp)
  colnames(ff) = c("actual", "forecast", "lower 80", "lower 95", "upper 80", "upper 95")
  cc__ = colnames(ff)
  colnames(ff_) <- cc__
  datatjes <- rbind(ff_, ff)
  colnames(datatjes) = c("actual", "forecast", "lower80", "lower95", "upper80", "upper95")
  datatjes <- data.frame(cbind(dates, datatjes))
  datatjes[,1] = dates
  colnames(datatjes)[2] = "cc"
  msfe = mean((cc_test- tmp[,1])^2)
  print(ggplot(data = datatjes, aes(x = dates, y = cc))+
         geom_line(col = "red") +
         geom_line(aes(y = forecast), col = "blue") +
         geom_line(aes(y = cc), col = "red") +
         geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.25) +
         ggtitle(paste("Consumer Confidence. Forecasted 6-step ahaead with \n Arima(", paste(sets[set,], collapse = ", "), ")", sep = "")) +
         xlab(paste("Dates", "     ", "E[MSFE]: ", as.character(round(msfe, 3)))) + 
         ylab("values"))
  if(set == 1){
    forecast.1step = tmp[,1]
  } else {
    forecast.1step = cbind(forecast.1step, tmp[,1])
  }
}
dev.off()

colnames(forecast.1step) = apply(format(sets), 1, paste, collapse=" ")
f = colnames(forecast.1step)
forecast.1step = cbind(cc_test, forecast.1step)
colnames(forecast.1step) = c("Actual", f)
write.csv(forecast.1step, "6 step ahead forecasts.csv")
write.csv(forecast.1step, "Forecasts Univariate.csv")









