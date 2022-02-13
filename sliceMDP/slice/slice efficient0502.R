#Dirichlet process prior
library(mvtnorm)
library(ggplot2)
set.seed(111)
n <- 100
p <- 2
nc0<-3
w0<- c(1/3,1/3,1/3)
cp0 <- sample(1:nc0,prob=w0,size=n,replace=TRUE)#component
mus0 <- cbind(c(-8,0,8),c(-8,0,8))
sds0 <- array(dim=c(p,p,nc0))
sds0[,,1]<-diag(p)
sds0[,,2]<-matrix(c(1,0.3,0.3,1),nrow=p)
sds0[,,3]<-matrix(c(1,0.8,0.8,1),nrow=p)
#y <- rnorm(n,mean=mus0[cp0],sd=sqrt(sds0[cp0]))
y <- matrix(ncol=p,nrow=n)
for (i in 1:n){
  y[i,] <- rmvnorm(1, mean = mus0[cp0[i],] , sigma = sds0[,,cp0[i]])
}
ggplot(as.data.frame(y), aes(x = y[,1], y = y[,2])) + geom_point(aes(col=cp0), size=3) + theme(legend.position="None")
#B. initial value####################################################
#prior for theta##############################B
mu0 <- rep(0,p)
k0 <- 0.1
v0 <- 2
T0 <- diag(p)
#prior for c
c<-1
#initial values
Nc <- n#######################################D
Ncy <- NULL
Nc_all <- NULL
delta <- c(1:n)#indicator variable############E
s <- exp(-delta)#deterministic decreasing sequence
theta <-  array(1, dim=c((p+1),p,n))#initial for theta###B
delta_all <- NULL
m_all <- NULL
v <- NULL
w <- NULL
u <- NULL#latent variable 1###################
#big loop##################################################################
I <- 1#total number of iterations
#big loop##################################################################
IN <- 1000
repeat{
  
  #1. updata oc,m
  oc <- vector("list",Nc)
  m <- NULL
  for (j in 1:Nc)
  {
    oc[[j]] <- which(delta==j)
    m[j] <- length(oc[[j]])
  }
  m_all <- c(m_all,length(which(m!=0)))#record all non-empty number of groups
  #2. update theta
  for (j in 1:Nc)
  {
    if (m[j]!=0){
      if ((m[j]==1)){
        ksi <- y[oc[[j]],]
        d <- matrix(0,nrow=p,ncol=p)
      }else{
        ksi <- colMeans(y[oc[[j]],])
        d <- matrix(0,nrow=p,ncol=p)
        for (i in 1:m[j]){
          d <- d+(y[oc[[j]][i],]-ksi)%*%t((y[oc[[j]][i],]-ksi))
        }
      }
      lambda1 <- v0+m[j]
      lambda2 <- solve(T0)+d+(k0*m[j]*(ksi-mu0)%*%t(ksi-mu0))/(k0+m[j])#inverse matrix
      mu1 <- (k0*mu0+m[j]*ksi)/(k0+m[j])
      mu2 <- k0+m[j]
      lambda <- rWishart(1,df=lambda1,Sigma=solve(lambda2))[,,1]#matrix
      mu <- rmvnorm(1,mean = mu1, sigma = solve(mu2*lambda))#vector
    } else {
      lambda <- rWishart(1,df=v0,Sigma=T0)[,,1]
      mu <- rmvnorm(1,mean = mu0, sigma = solve(k0*lambda))
    }
    theta[,,j]<-rbind(mu,lambda)#results from last iteration
  }
  #3. update v
  for (j in 1:(Nc-1)){
    v[j] <- rbeta(1,(1+m[j]),(c+cumsum(m[(j+1):Nc])))
  }
  v[Nc] <- rbeta(1,(1+m[Nc]),c)
  
  #4. update u
  for (ia in 1:n){
    u[ia] <- runif(1,0,s[delta[ia]])
    Ncy[ia] <- max(which(u[ia]<s))
  }
  Nc <- max(Ncy)#new number of groups
  Nc_all <- c(Nc_all, Nc)#record all max number of groups
  
  #5. update w  
  w[1]<-v[1]
  for (id in (2:Nc)){
    w[id]<-v[id]*prod(1-v[1:(id-1)])
  }
  
  #6. update delta
  prb<-matrix(nrow=n,ncol=Nc)
  for (ie in 1:n){
    for (je in (1:Nc)){
      if (s[je]>u[ie]){
        prb[ie,je]<-w[je]*dmvnorm(y[ie,],mean=theta[1,,je],sigma=solve(theta[-1,,je]))/s[je]
        }else{
        prb[ie,je]<-0
      }
    }
  }
  delta <- max.col(prb)
  delta_all <- rbind(delta_all, delta)#record all group allocation

  I<-I+1
  if (I>IN)
    break
}

m_all[400:500]

cluster<-m_all[-1:-500]
average<-NULL
for (i in 1:length(cluster)){
  average[i]<-mean(cluster[1:i])
}
average[length(cluster)]

