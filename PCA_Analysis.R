
###################################################################################
# CONSUMER CONFIDENCE USING GOOGLE TRENDS
# Author: Tom Kennes           
# Date: 03/12/2017           
###################################################################################

# Clear workspace memory
rm(list=ls())
cat("\014")

# Set working directory and load data
setwd("####################")
rawX <- read.csv('GoogleQueries.csv', header = TRUE, sep = ';')
rawy <- read.csv('CBS_cc.csv', header = TRUE, sep = ';')
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
    dates = as.Date(dates, tz = "GMT")
    return(dates)
  } else {
    colnames(dats) = c("years", "months")
    return(dats)
  }
}

# Set sizes
n <- nrow(rawX)
p <- ncol(rawX)
# First two columns of X are irrelevant (variable time)
rawX <- rawX[ , 3:p]
rawy <- as.matrix(rawy)
p <- p - 2

# Testing for unit roots
library(tseries)
adf.test(rawy[!is.na(rawy)], alternative = "stationary", k = 0)
adf <- vector("numeric", p)
for(i in 1:p)
{
  adf[i] <- adf.test(rawX[,i])$p.value
}

# Difference the explanatory variables if nonstationary
X <- matrix(NA, nrow = 118, ncol = 88)
for(i in 1:p) {
  if(adf[i] < 0.05) {
    X[,i] <- rawX[2:119,i]
  }
  else if (adf[i] >= 0.05) {
    X[,i] <- diff(as.matrix(rawX[,i]))
  }
}

# Difference the consumer confidence index if nonstationary
if(adf.test(rawy)$p.value > 0.01)
{
  y <- diff(rawy)
}
n <- nrow(X)

# Split sample into training (estimation) and testing (validation)
Xest <- X[1:95,]
yest <- y[1:95,]
Xval <- X[96:118,]
yval <- y[96:118,]

# Standardize the variables
colSd <- function (x, na.rm=FALSE) apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)
Xmeanest <- colMeans(Xest)
Xsdest <- colSd(Xest)
ymeanest <- mean(yest)
ysdest <- sd(yest)
Xest <- scale(Xest, TRUE, TRUE)
yest <- scale(yest, TRUE, TRUE)
Xval <- scale(Xval,center = Xmeanest, scale = Xsdest)
yval <- scale(yval,center = ymeanest, scale = ysdest)


###################################################################################
# DETERMINING THE NUMBER OF FACTORS     
###################################################################################


# Determining the number of factors: Information Criteria
library(POET)
POETKhat(Xest)

# Determining the number of factors: Eigenvalues
library(nFactors)
ev <- eigen(cor(Xest))                                                          # get eigenvalues
dist <- parallel(subject = nrow(Xest),var = ncol(Xest), rep = 100, cent = .05)  # distribution of the eigenvalues
nF <- nScree(x = ev$values, aparallel = dist$eigen$qevpea)                      # returns analysis of number of factors
plotnScree(nF)


###################################################################################
# FACTOR ANALYSIS AND PRINCIPAL COMPONENTS         
###################################################################################


# Principal Component Analysis
library(factoextra)
pca1 <- prcomp(Xest, retx=TRUE, center=TRUE, scale=TRUE)

# Percentage of Variance Explained
pr_var <- ((pca1$sdev)^2)                                                       # Compute the variance of the principal components
pr_var[1:20]                                                                    # Check the variance of the first 10 principal components
prop_varex <- pr_var/sum(pr_var)   # Compute proportion of variance explained
library(R6)
fviz_eig(pca1, ncp = 25, xlab = "Principal Component", ylab = "Percentage of Variance Explained", main ="")                               
total_var <- sum(pr_var[1:7])                                                   # Compute total variance explained by principal components                 
number_factors <- 10

# Set up the Training model
lag_yest <- c(NA, yest[1:94])                                                   # Create lag of CCI
train.data <- data.frame(CCI = yest, CCI_lag = lag_yest, pca1$x)                # Create training dataset
train.data <- train.data[,1:(number_factors+2)]                                 # Take the first r principal components
linear.model <- lm(CCI ~ ., data = train.data)
summary(linear.model)


###################################################################################
# PREDICTIONS 23 STEP AHEAD 
###################################################################################


n = nrow(Xval)
predictions <- rep(0, n)
factor <- (Xval[1,] %*% pca1$rotation)[1:(number_factors+2)]
predictions[1] = linear.model$coefficients[1] +
  linear.model$coefficients[2]*train.data[nrow(train.data),1] +
  linear.model$coefficients[3:length(linear.model$coefficients)] %*% factor[3:(number_factors+2)]
for(i in 2:n){
  factor <- (Xval[i,] %*% pca1$rotation)[1:(number_factors+2)]
  predictions[i] = linear.model$coefficients[1] +
    linear.model$coefficients[2]*predictions[i - 1]+
    linear.model$coefficients[3:length(linear.model$coefficients)] %*% factor[3:(number_factors+2)]
}

# Plotting the predictions
dat1 = c(2008, 1)
dat2 = c(2017, 11)
dates = dates_creator(dat1, dat2)
dates <- dates[2:(length(dates) - 1)]
E <- c(linear.model$fitted.values, predictions)
datatjes <- data.frame(dates, y[2:length(y)], E)
require(ggplot2)
ggplot(data = datatjes, aes(x = dates, y = y.2.length.y..)) + 
  geom_line(col = "red") + 
  geom_line(aes(y = E), col = "blue") +
  ggtitle("") + 
  ylab("Values") + 
  xlab("Time") +
  geom_vline(aes(xintercept = as.numeric(datatjes[95,1])), linetype = 4)


# Plotting the predictions
dat1 = c(2008, 1)
dat2 = c(2017, 11)
dates = dates_creator(dat1, dat2)
dates <- dates[2:(length(dates) - 1)]
E <- c(linear.model$fitted.values, predictions)
y <- c(yest, yval)
datatjes <- data.frame(dates, y[2:length(y)], E)
require(ggplot2)
ggplot(data = datatjes, aes(x = dates, y = y.2.length.y..)) + 
  geom_line(col = "red") + 
  geom_line(aes(y = E), col = "blue") +
  ggtitle("") + 
  ylab("Values") + 
  xlab("Time") +
  geom_vline(aes(xintercept = as.numeric(datatjes[95,1])), linetype = 4)


# Compute Mean Squared Forecast Error
library(Metrics)
mse(yval, predictions)


###################################################################################
# PREDICTIONS 6 STEP AHEAD 
###################################################################################


# Rolling Time Window (6-step ahead)
predictions.ahead <- function(Xest, yest, ahead, number_of_factors){
  number_factors = number_of_factors
  pca1 <- prcomp(Xest, retx=TRUE, center=TRUE, scale=TRUE)
  lag_yest <- c(NA, yest[1:(length(yest)-1)])                                                   # Create lag of CCI
  train.data <- data.frame(CCI = yest, CCI_lag = lag_yest, pca1$x)                # Create training dataset
  train.data <- train.data[,1:(number_factors+2)]                                 # Take the first r principal components
  linear.model <- lm(CCI ~ ., data = train.data)
  n = ahead
  predictions <- rep(0, n)
  factor <- (Xval[1,] %*% pca1$rotation)[1:(number_factors+2)]
  predictions[1] = linear.model$coefficients[1] +
    linear.model$coefficients[2]*train.data[nrow(train.data),1] +
    linear.model$coefficients[3:length(linear.model$coefficients)] %*% factor[3:(number_factors+2)]
  for(i in 2:n){
    factor <- (Xval[i,] %*% pca1$rotation)[1:(number_factors+2)]
    predictions[i] = linear.model$coefficients[1] +
      linear.model$coefficients[2]*predictions[i - 1]+
      linear.model$coefficients[3:length(linear.model$coefficients)] %*% factor[3:(number_factors+2)]
  }
  return(predictions)
}  

# 6 Step ahead predictions
to.predict.ahead = 6
forecast <- rep(0, n)
for(i in 0:(nrow(Xval)-to.predict.ahead)){
  if(i == 0){
    forecast[1:to.predict.ahead] <- predictions.ahead(Xest, yest, ahead=to.predict.ahead, number_of_factors = 7) 
  } else {
    ynew <- c(yest, yval[1:i, ])
    xnew <- rbind(Xest, Xval[1:i, ])
    tmp  <- predictions.ahead(xnew, ynew, ahead=to.predict.ahead, number_of_factors = 7) 
    forecast[i + to.predict.ahead] <- tmp[to.predict.ahead]
  }
}

# Plotting the predictions
dat1 = c(2008, 1)
dat2 = c(2017, 11)
dates = dates_creator(dat1, dat2)
dates <- dates[2:(length(dates) - 1)]
E <- c(linear.model$fitted.values, forecast)
y <- c(yest, yval)
datatjes <- data.frame(dates, y[2:length(y)], E)
require(ggplot2)
ggplot(data = datatjes, aes(x = dates, y = y.2.length.y..)) + 
  geom_line(col = "red") + 
  geom_line(aes(y = E), col = "blue") +
  ggtitle("") + 
  ylab("Values") + 
  xlab("Time") +
  geom_vline(aes(xintercept = as.numeric(datatjes[95,1])), linetype = 4)

# Compute Mean Squared Forecast Error
library(Metrics)
mse(yval, predictions)



