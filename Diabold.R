###################################################################################
# CONSUMER CONFIDENCE USING GOOGLE TRENDS
# Author: Tom Kennes            
# Date: 03/12/2017           
###################################################################################

require("forecast")
setwd("#####################################")

factor <- read.csv("Forecasts Factor Model.csv", header = T, sep = ";")
c <- colnames(factor)[2:ncol(factor)]
lasso <- read.csv("Lasso Forecasts.csv", header =T, sep = ";", dec = ",")
uni <- read.csv("Forecasts Univariate.csv", sep = ",")
lasso <- lasso[,2:ncol(lasso)]
colnames(lasso) <- c("ADALasso", "Post Lasso", "Post Double Selection")

factor <- factor[, 2:ncol(factor)]


uni <- uni[,which(colnames(uni) %in% c("Actual", "X5.0.5", "X1.0.0"))]
colnames(uni) =  c("Actual", "ARIMA5..0..5", "ARIMA1..0..0")
uni = uni[,c(1,3,2)]
observations = uni$Actual
forecasts <- cbind(uni, factor, lasso)
mse = matrix(0, ncol = 1, nrow = ncol(forecasts))

#plot(lasso$Observed, col = "red", type = "l", ylim = c(max(uni$Actual), min(uni$Actual)))
#points(seq(1, length(uni$Actual), by = 1), uni$Actual, col = "blue", type = "l")


for(i in 1:ncol(forecasts)){
  mse[i,1] = mean((observations - forecasts[,i])^2)
}
row.names(mse) = colnames(forecasts)
mse = t(t(mse[2:nrow(mse),]))
mse = t(t(mse[order(mse[,1]),]))
mse = t(t(mse[2:nrow(mse),]))
mse = t(t(mse[order(mse[,1]),]))
colnames(mse) <- c("MSE")
write.csv(mse, "mse.csv")

forecasts <- forecasts[,2:ncol(forecasts)]
diabold_mariano <- matrix(0, nrow = ncol(forecasts), ncol = ncol(forecasts))
for(i in 1:ncol(forecasts)){
  for(j in 1:ncol(forecasts)){
    if(i == j){
      diabold_mariano[i,j] = 1
    } else {
      diabold_mariano[i,j] = dm.test(forecasts[,i], forecasts[,j], alternative = "two.sided", h = 22)$p.value
    }
  }
}



colnames(diabold_mariano) <- colnames(forecasts)
row.names(diabold_mariano) <- colnames(forecasts)
write.csv(diabold_mariano, "diabold.csv")
diabold_mariano
