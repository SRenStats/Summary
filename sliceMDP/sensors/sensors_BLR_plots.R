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
#3.prior information###########################################
mu<-colMeans(beta_training)#mean vector for beta
s<-NULL
for (i in 1:p)
{
  s[i]<-crossprod(beta_training[,i]-mu[i])/(nrow(beta_training)-1)
}
sigma<-diag(s)#diagnal variance matrix for beta
lambda<-mean(1/delta2_training)#for y
#4. update beta using data(x,y)#######################################
number_test<-dim(data_test)[3]
beta_test<-NULL
delta2_test<-NULL
for (i in 1:number_test){
  x<-cbind(1,data_test[,-1,i])
  y<-data_test[,1,i]
  a<-solve(lambda*t(x)%*%x+solve(sigma))%*%(lambda*t(x)%*%y+solve(sigma)%*%mu)
  b<-crossprod(y-x%*%a)/(length(y)-p)
  beta_test<-rbind(beta_test,t(a))
  delta2_test<-c(delta2_test,b)
}

###5.plots for delta2
###1-red dots for all ols
###2-green dots for training data ols
###3-blue dots for baysian linear regression
#5.1 scatter plots
p1<-cbind(1,c(1:length(delta2_all)),delta2_all)
p2<-cbind(2,c(1:length(delta2_training)),delta2_training)
p3<-cbind(3,c((length(delta2_training)+1):length(delta2_all)),delta2_test)
plotdata<-data.frame(
  x=c(p1[,2],p2[,2],p3[,2]),
  y=c(p1[,3],p2[,3],p3[,3]),
  c=factor(c(p1[,1],p2[,1],p3[,1]))
)
par(mfrow=c(1,1))
ggplot(plotdata, aes(x =x, y = y,colour=c)) + geom_point() + labs(title=name,x="sensors",y="delta2") + theme(plot.title = element_text(hjust = 0.5))

#5.2 density plots
par(mfrow=c(1,1))
plot(density(delta2_training),col="green",main=name)
lines(density(delta2_all),col="red",main=name)
lines(density(delta2_test),col="blue")
lines(density(delta2_all[(length(delta2_training)+1):length(delta2_all)]),col="black",main=name)

#5.3 line plots
par(mfrow=c(2,1))
plot(c(1:length(delta2_all)),delta2_all,col="red",main=name,xlab="sensors",ylab="delta2")
#lines(c(1:length(delta2_training)),delta2_training,col="green")
lines(c((length(delta2_training)+1):length(delta2_all)),delta2_test,col="blue")
plot(c(1:length(delta2_all)),(delta2_all),col="red",main=name,xlab="sensors",ylab="delta2")
lines(c((length(delta2_training)+1):length(delta2_all)),delta2_all[(length(delta2_training)+1):length(delta2_all)],col="black")

#5.4 plots for beta
par(mfrow=c(1,5))
for (j in 1:5){
plot(density((beta_training[,j])),col="black",main=c(name,j))
}
