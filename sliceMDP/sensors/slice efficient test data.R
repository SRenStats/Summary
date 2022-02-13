#Dirichlet process prior
library(mvtnorm)
library(ggplot2)
set.seed(111)
n <- 100
x <- cbind(1,rmvnorm(n, mean = c(2,3,4,5) , sigma = diag(4)))#fixed x
p <- 5
nc0 <- 3#assume beta0 comes from a mixture of 3
w0 <- c(0.1,0.3,0.6)
cp0 <- sample(1:nc0,prob=w0,size=n,replace=TRUE)#component
mus0 <- matrix(nrow=nc0,ncol=p)#3*5 matrix for 3 clusters and 5 dimensions
mus0[1,] <- c(2,3,4,5,6)#mean for cluster 1
mus0[2,] <- c(-2,-3,-4,-5,-6)
mus0[3,] <- c(2,-3,4,-5,6)
sds0 <- array(dim=c(p,p,nc0))#5*5*3 array for 5 dimension covariance matrix and 3 clusters
sds0[,,1]<-diag(p)#5*5 matrix
sds0[,,2]<-diag(p)
sds0[,,3]<-diag(p)

beta0 <- matrix(nrow=n,ncol=p)
for (i in 1:n){
  beta0[i,] <- rmvnorm(1, mean = mus0[cp0[i],] , sigma = sds0[,,cp0[i]])
}


#ep0 <- rnorm(n)
#y <- x%*%beta0 + ep0#fixed x, DPM beta



