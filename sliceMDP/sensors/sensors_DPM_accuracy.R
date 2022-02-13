library(ggplot2)
library(R.matlab)
data<-readMat("D:/phd/sensors/data/T06.mat")$Data
data_rdc<-readMat("D:/phd/sensors/data/T06_Reduced.mat")$Data2
number_sensors<-dim(data)[3]
p<-5

#1.estimate beta using all data#################################
beta_all<-NULL
delta2_all<-NULL
for (i in 1:number_sensors){
  x<-cbind(1,data[,-1,i])
  y<-data[,1,i]
  a<-solve(t(x)%*%x)%*%t(x)%*%y#OLS estimate
  b<-crossprod(y-x%*%a)/(length(y)-p)
  beta_all<-rbind(beta_all,t(a))
  delta2_all<-c(delta2_all,b)
}
#2.estimate beta using training data#################################
#number_training<-floor(number_sensors/2)#integer
#data_training<-data[,,1:number_training]
number_training<-240#integer
data_training<-data[,,1:number_training]
data_test<-data_rdc[,,(number_training+1):number_sensors]
data_accu<-data[,,(number_training+1):number_sensors]
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
#3.DPM as prior information###########################################
#In the former method, we use sample mean and sample variance of beta as mean and variance of prior.
#Here, we asumme beta from DPM instead. beta~w*N(mu,sigma)

unique(delta)
k_s <- length(unique(delta))
mu_s <- matrix(nrow=p,ncol=k_s)
sigma_s <- array(dim=c(p,p,k_s))
oc_s <- vector("list",k_s)
m_s <- NULL

for (j in 1:k_s)
{
  oc_s[[j]] <- which(delta==unique(delta)[j])
  m_s[j] <- length(oc_s[[j]])
#}
  beta_training2 <- beta_training[oc_s[[j]],]
  mu<-colMeans(beta_training2)#mean vector for beta
  s<-NULL
  for (i in 1:p)
  {
    s[i]<-crossprod(beta_training2[,i]-mu[i])/(nrow(beta_training2)-1)
  }
  sigma<-diag(s)#diagnal variance matrix for beta
  mu_s[,j] <- mu
  sigma_s[,,j] <- sigma
}
w_s <- m_s/number_training
#4. update beta using training data and lambda, mu, sigma####
number_test<-dim(data_test)[3]
beta_test<-NULL
delta2_test<-NULL
delta2 <- NULL#new allocation
lambda <- 1/mean(delta2_training)
for (i in 1:number_test){
  x<-cbind(1,data_test[,-1,i])
  y<-data_test[,1,i]
  prb_k <- NULL#probability of each sensor for each cluster
  for (j in 1:k_s){
    prb_k[j] <- w_s[j]*dmvnorm(y,mean=x%*%mu_s[,j],sigma=mean(delta2_training)*diag(length(y)))    
  }
  delta2[i] <- which.max(prb_k)#new allocation
  mu <- mu_s[,delta2[i]]
  sigma1 <- sigma_s[,,delta2[i]]
  a<-solve(lambda*t(x)%*%x+solve(sigma1))%*%(lambda*t(x)%*%y+solve(sigma1)%*%mu)
  b<-crossprod(y-x%*%a)/(length(y)-p)
  beta_test<-rbind(beta_test,t(a))
  delta2_test<-c(delta2_test,b)
}

#5.estimate accuracy#####
c <- NULL
for (i in 1:number_test){
  x<-cbind(1,data_accu[,-1,i])
  y<-data_accu[,1,i]
  a <- abs(y-x%*%beta_test[i,]) 
  b <- c(mean(a[1:6]),mean(a[7:12]),mean(a[13:18]),mean(a[19:26]),mean(a[27:32]))
  c <- c(c,length(which(b>0.003)))#number of not passing
}
length(which(c==0))/number_test#pass rate, the higher the better

