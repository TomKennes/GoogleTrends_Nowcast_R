#Installing necessary packages
install.packages("glmnet")
install.packages("tseries")

###################################################################################
# CONSUMER CONFIDENCE USING GOOGLE TRENDS
# Author: Tom Kennes            
# Date: 09/12/2017 
# Results may be slightly different due to randomness in the cross validation
##################################################################################

## READ AND CLEAN DATA

#Clear workspace memory
rm(list=ls())
set.seed = 1234

#Set working directory and load data
setwd("####################")
rawX <- read.csv2("GoogleQueries01.csv", check.names = FALSE)
rawy <- read.csv2("CBS_cc.csv", header = FALSE)

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
y <- y[2:n]
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
Xfull <- rbind(Xest,Xval)
yfull <- rbind(yest,yval)
  
## ADAPTED LASSO

library("glmnet")
#Do ridge regression to obtain initial estimates for beta
ridge.lambda <- vector("numeric", 100)
for(i in 1:100)
{
  cv.ridge <- cv.glmnet(Xest[,2:p], yest, alpha = 0, nfolds = 10)
  ridge.lambda[i] <- cv.ridge$lambda.min
}
ridge <- glmnet(Xest[,2:p],yest, alpha = 0, lambda = mean(ridge.lambda))
ridge.beta <- as.matrix(ridge$beta)
#create weights vector
w <- as.matrix(c(0,1/abs(ridge.beta)))
#determine the optimal value for lambda
adplasso.lambda <- vector("numeric", 100)
for(i in 1:100)
{
  cv.adplasso <- cv.glmnet(x = Xest,y = yest, family="gaussian", penalty.factor = w, nfolds = 10)
  adplasso.lambda[i] <- cv.adplasso$lambda.min
}
#build lasso model according to this lambda
adplasso <- glmnet(Xest, yest, family = "gaussian", penalty.factor = w, alpha = 1, lambda = 0.55*mean(adplasso.lambda))
adplasso.beta <- as.matrix(adplasso$beta)
#do prediction with a rolling window, 23 steps ahead
pred.y <- vector("numeric",nrow(Xval))
pred.y[1] <- Xval[1,]%*%adplasso.beta
for(i in 2:nrow(Xval))
{
  pred.y[i] <- cbind(pred.y[i-1],t(as.vector(Xval[i,2:ncol(Xval)])))%*%adplasso.beta
}
#do prediction with a six step ahead forecast rolling window
pred.y.6step <- vector("numeric",nrow(Xval))
for(i in -6:(nrow(Xval)-6))
{
  #Estimate new model with additional data point, to determine T + i + 6
  adplasso <- glmnet(Xfull[1:(nrow(Xest)+i),], yfull[1:(nrow(yest)+i)], family = "gaussian", penalty.factor = w, alpha = 1, lambda = 0.55*mean(adplasso.lambda))
  adplasso.beta <- as.matrix(adplasso$beta)
  #Do six step ahead forecast with this model
  sixstep <- vector("numeric",6)
  sixstep[1] <- Xfull[(nrow(Xest)+i+1),]%*%adplasso.beta
  for(j in (i+2):(i+6))
  {
    sixstep[j-i] <- cbind(sixstep[j-i-1],t(as.vector(Xfull[nrow(Xest)+j,2:ncol(Xval)])))%*%adplasso.beta
  }
  #Store 6 step ahead forecast
  pred.y.6step[i+6] <- sixstep[6]
}
  
## POST-LASSO

#create matrix with variables corresponding to non-zero coefficients
X.post.est <- matrix(0,nrow(Xest),0)
X.post.val <- matrix(0,nrow(Xval),0)
coef.name <- matrix(0,0,1)
for(j in 1:p)
{
  if(adplasso.beta[j] != 0)
  {
    X.post.est <- cbind(X.post.est,Xest[,j])
    X.post.val <- cbind(X.post.val,Xval[,j])
    coef.name <- rbind(coef.name, rownames(adplasso.beta)[j])
  }
}
X.full.post <- rbind(X.post.est,X.post.est)
#run OLS on relevant variable matrix
b.post <- solve(crossprod(X.post.est))%*%crossprod(X.post.est,yest)
rownames(b.post) <- coef.name
#predict using these values for beta, on a rolling window, 23 steps ahead
pred.y.post <- vector("numeric",nrow(Xval))
pred.y.post[1] <- X.post.val[1,]%*%b.post
for(i in 2:nrow(X.post.val))
{
  pred.y.post[i] <- cbind(pred.y.post[i-1],t(as.vector(X.post.val[i,2:ncol(X.post.val)])))%*%b.post
}
#do prediction with a six step ahead forecast rolling window
pred.y.6step.post <- vector("numeric",nrow(Xval))
pred.y.6step.post[1] <- pred.y.post[6]
for(i in -6:(nrow(Xval) - 6))
{
  #run OLS on relevant variable matrix, up to information in T+i
  b.post <- solve(crossprod(X.full.post[1:(nrow(Xest)+i),]))%*%crossprod(X.full.post[1:(nrow(Xest)+i),],yfull[1:(nrow(yest)+i)])
  #Do six step ahead forecast with this model
  sixstep <- vector("numeric",6)
  sixstep[1] <- X.full.post[(nrow(Xest)+i+1),]%*%b.post
  for(j in (i+2):(i+6))
  {
    sixstep[j-i] <- cbind(sixstep[j-i-1],t(as.vector(X.full.post[nrow(Xest)+j,2:ncol(X.post.val)])))%*%b.post
  }
  #Store 6 step ahead forecast
  pred.y.6step.post[i+6] <- sixstep[6]
}

## POST DOUBLE SELECTION

#Do post-double selection selection, to do inference
pds.inf <- matrix(0,nrow(b.post),4)
rownames(pds.inf) <- coef.name
colnames(pds.inf) <- c("beta","standard error","t-stat","p-value")
for(i in 1:nrow(adplasso.beta))
{
  if(adplasso.beta[i] != 0)
  {
    Xpd <- Xest[,-i]
    Wpd <- as.matrix(Xest[,i])
    penaltypd <- w[-i]
    #regress W on X and Y on X
    delta <- as.matrix(glmnet(Xpd, yest, family = "gaussian", penalty.factor = penaltypd, alpha = 1, lambda = 0.55*mean(adplasso.lambda))$beta)
    gamma <- as.matrix(glmnet(Xpd, Wpd, family = "gaussian", penalty.factor = penaltypd, alpha = 1, lambda = 0.55*mean(adplasso.lambda))$beta)
    #Create set of at least one selected regressor
    XS.hat <- matrix(0,nrow(Xest),0)
    for(j in 1:ncol(Xpd))
    {
      if((delta[j] != 0 )||(gamma[j] != 0))
      {
        XS.hat <- cbind(XS.hat,Xpd[,j])
      }
    }
    XS.hat <- cbind(XS.hat, Wpd)
    #Regress Y on XS.hat, using standard OLS
    beta.pds <- solve(crossprod(XS.hat))%*%crossprod(XS.hat, yest)
    sigma.hat <- crossprod(yest - XS.hat%*%beta.pds)/(nrow(XS.hat)-ncol(XS.hat))
    se.beta <- sqrt(solve(crossprod(XS.hat))[[ncol(XS.hat),ncol(XS.hat)]]*sigma.hat)
    if(beta.pds[nrow(beta.pds)] > 0)
    {
      pval <- 2*pnorm(beta.pds[nrow(beta.pds)]/se.beta,lower.tail = FALSE)
    }
    else
    {
      pval <- 2*pnorm(beta.pds[nrow(beta.pds)]/se.beta,lower.tail = TRUE)
    }
    pds.inf[row.names(adplasso.beta)[i],1] <- beta.pds[nrow(beta.pds)]
    pds.inf[row.names(adplasso.beta)[i],2] <- se.beta
    pds.inf[row.names(adplasso.beta)[i],3] <- beta.pds[nrow(beta.pds)]/se.beta
    pds.inf[row.names(adplasso.beta)[i],4] <- pval
  }
}
#Predict using pds betas
pred.y.pds <- vector("numeric",nrow(Xval))
pred.y.pds[1] <- X.post.val[1,]%*%pds.inf[,1]
for(i in 2:nrow(Xval))
{
  pred.y.pds[i] <- cbind(pred.y.pds[i-1],t(as.vector(X.post.val[i,2:ncol(X.post.val)])))%*%pds.inf[,1]
}

## PLOT THE RESULTS

#Using 23 step ahead
plot(yval, type = "l", col = "black", xlab = "Time", ylab = "Difference in Consumer Confidence", main = "23-step ahead forecast")
lines(pred.y.post, type = "l", col = "red")
lines(pred.y, type = "l", col = "blue")
lines(pred.y.pds, type = "l", col = "green")
legend(15,1.4,c("True data","Post-Lasso","Adapted Lasso","Post Double Selection"),col=c("black","red","blue","green"),lty=1)
#Using rolling window with 6 step ahead
plot(yval, type = "l", col = "black", xlab = "Time", ylab = "Difference in Consumer Confidence", main = "6-step ahead forecast with rolling window")
lines(pred.y.6step.post, type = "l", col = "red")
lines(pred.y.6step, type = "l", col = "blue")
legend(15,1.4,c("True data","Post-Lasso","Adapted Lasso"),col=c("black","red","blue"),lty=1)

#Calculate MSE
MSE.adplasso <- (1/nrow(Xval))*sum((yval-pred.y)^2)
MSE.post <- (1/nrow(Xval))*sum((yval-pred.y.post)^2)
MSE.pds <- (1/nrow(Xval))*sum((yval-pred.y.pds)^2)
MSE <- matrix(c(MSE.adplasso, MSE.post, MSE.pds))
rownames(MSE) <- c("Adaptive Lasso","Post-Lasso","Post-Double Selection")

#Write the forecasts to a CSV
write.csv(cbind(yval,pred.y.6step,pred.y.6step.post), file = "Lasso Forecasts.csv", sep = ";")
