library(ggplot2)
library(R.matlab)
big<-readMat("D:/sensors/data/T25.mat")$Data
big_rdc<-readMat("D:/sensors/data/T25_Reduced.mat")$Data2
medium<-readMat("D:/sensors/data/T15.mat")$Data
medium_rdc<-readMat("D:/sensors/data/T15_Reduced.mat")$Data2
small<-readMat("D:/sensors/data/T06.mat")$Data
small_rdc<-readMat("D:/sensors/data/T06_Reduced.mat")$Data2
data<-small#
data_rdc<-small_rdc
name<-"small"#title for plots
number_sensors<-dim(data)[3]
p<-5

#1.estimate beta using all data#################################
beta_all<-NULL
delta2_all<-NULL
for (i in 1:number_sensors){
  x<-cbind(1,data[,-1,i])
  y<-data[,1,i]
  a<-solve(t(x)%*%x)%*%t(x)%*%y
  b<-crossprod(y-x%*%a)/(length(y)-p)
  beta_all<-rbind(beta_all,t(a))
  delta2_all<-c(delta2_all,b)
}
#2.estimate beta using training data#################################
number_training<-240#integer
data_training<-data[,,1:number_training]
data_test<-data_rdc[,,(number_training+1):number_sensors]
beta_training<-NULL
delta2_training<-NULL
for (i in 1:number_training){
  x<-cbind(1,data_training[,-1,i])
  y<-data_training[,1,i]
  a<-solve(t(x)%*%x)%*%t(x)%*%y
  b<-crossprod(y-x%*%a)/(length(y)-p)
  beta_training<-rbind(beta_training,t(a))
  delta2_training<-c(delta2_training,b)
}
#3 plots for beta########
j<-0
j <- j+1
hist(beta_all[,j]*10^4,breaks=100,col="blue1",main=c(name,j))
sort(beta_all[,j]*10^4)
plot(density((beta_all[,j])),col="black",main=c(name,j))
