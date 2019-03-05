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



x= c(2008,1)
y = c(2017, 11)
dates = dates_creator(x, y)
cc_training = cc[1:95]
cc_test = cc[96:length(cc)]

require(ggplot2)
ggplot( data = data, aes( dates, cc )) + 
  geom_line(size = 2, color = "#0072B2") + 
  ggtitle("Consumer Confidence") + 
  xlab("Time") + 
  ylab("Value")
