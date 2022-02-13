#1.generate data from a mixture model#########################################
library(MASS)
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
y <- matrix(ncol=p,nrow=n)
for (i in 1:n){
  y[i,] <- rmvnorm(1, mean = mus0[cp0[i],] , sigma = sds0[,,cp0[i]])
}
ggplot(as.data.frame(y), aes(x = y[,1], y = y[,2])) + geom_point(aes(col=cp0), size=3) + theme(legend.position="None")
#2. initial value####################################################
#prior for theta##############################B
mu0 <- rep(0,p)
k0 <- 0.1
v0 <- 2
T0 <- diag(p)
#prior for c
a<-0.1
b<-100
c<-0.5
#initial values
Nc <- n#######################################D
Ncy <- NULL
Nc_all <- NULL
delta <- c(1:n)#indicator variable############E
#theta <-  matrix(1, nrow=2, ncol=n)#initial for theta###B
theta <-  array(1, dim=c((p+1),p,n))#initial for theta###B
delta_all <- NULL
m_all <- NULL
v <- rbeta(n,1,c)# 1*c/((1+c+1)*(1+c)^2)
w <- NULL
w[1] <- v[1]
for (i in (2:n)){
  w[i] <- v[i]*prod(1-v[1:(i-1)])}
u <- NULL#latent variable 1###################A
#big loop#################################################################
I<-1
#big loop##################################################################
IN<-1000
repeat{ 
  #A.update u#######################################################
  for (ia in 1:n){
    u[ia] <- runif(1,0,w[delta[ia]])
  }
  #observations within each component
  oc<-vector("list",Nc)
  m<-NULL
  for (j in 1:Nc)
  {
    oc[[j]]<-which(delta==j)
    m[j] <-length(oc[[j]])
  }
  m_all <- c(m_all,length(which(m!=0)))#number of non-empty groups
  #B.update theta=(mu,lambda)###################################################
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
  #C.update v##########################################
  alpha<-NULL
  beta<-NULL
  for (j in 1:Nc)
  {
    if (j==1){
      alpha0 <- u[oc[[j]]]
    }else{
      alpha0 <- u[oc[[j]]]/prod(1-v[1:(j-1)])
    }
    if (length(alpha0)==0){
      alpha[j] <- 0
    } else {
      alpha[j] <- max(alpha0)
    }
    beta0<-NULL
    if (j<Nc)
    {
      for (j2 in (j+1):Nc)
      {
        beta0<-c(beta0,u[oc[[j2]]]*(1-v[j])/(v[j2]*prod(1-v[1:(j2-1)])))
      }
      if (length(beta0)==0){
        beta[j]<-1
      } else {
        beta[j]<-1-max(beta0) 
      }
    }else{
      beta[j]<-1
    }
    if (alpha[j]<beta[j])
    {  unif<-runif(1)
    v[j]<-1-((1-alpha[j])^c*(1-unif)+(1-beta[j])^c*unif)^(1/c)
    }
  }  
  #D.update w and the number of components Nc###########################################
  sum_w<-NULL
  w[1]<-v[1]
  sum_w[1]<-w[1]
  for (id in (2:n)){
    w[id]<-v[id]*prod(1-v[1:(id-1)])
    sum_w[id]<-sum(w[1:id])
  }
  for (id2 in 1:n){
    if(length(which(sum_w>1-u[id2]))!=0)
    { Ncy[id2] <- min(which(sum_w>1-u[id2]))###correct
    } else {
      Ncy[id2] <- n
    }
  }
  Nc <- max(Ncy)#new number of groups
  Nc_all<-c(Nc_all,Nc)
  #E.update delta######################################################################
  prb<-matrix(nrow=n,ncol=Nc)
  for (ie in 1:n){
    for (je in (1:Nc)){
      if (w[je]>u[ie]){
        prb[ie,je]<-dmvnorm(y[ie,],mean=theta[1,,je],sigma=solve(theta[-1,,je]))
      }else{
        prb[ie,je]<-0
      }
    }
  }
  delta<-max.col(prb)
  delta_all<-rbind(delta_all,delta)#group allocation
  #F.update c##########################################################################gibbs sampling
  k<-length(unique(delta))
  c_g<-NULL
  eta<-NULL
  eta[1]<-rbeta(1,c+1,n)
  nf<-5000
  for (iff in 2:nf){
    pi<-(a+k-1)/(n*b-n*log(eta[iff-1])+a+k-1)
    unif<-runif(1)
    if (unif<pi){
      c_g[iff]<-rgamma(1,shape=a+k,rate=b-log(eta[iff-1]))
    }else{
      c_g[iff]<-rgamma(1,shape=a+k-1,rate=b-log(eta[iff-1]))
    }
    eta[iff]<-rbeta(1,c_g[iff]+1,n)
  }
  c<-mean(c_g[(nf/2):nf])
  I<-I+1
  if (I>IN)
    break
}
#View(delta_all)
mean(log(eta))
m_all
sort(m[which(m!=0)])#number of components in each cluster
cluster<-m_all[-1:-IN/2]
average<-NULL
for (i in 1:length(cluster)){
  average[i]<-mean(cluster[1:i])
}
average[length(cluster)]
