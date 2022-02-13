#Dirichlet process prior
#A.generate data from a mixture model#########################################
set.seed(111)
n <- 50
nc0 <- 3
w0<- c(1/3,1/3,1/3)
cp0 <- sample(1:nc0,prob=w0,size=n,replace=TRUE)#component
mus0 <- c(-8,0,8)
sds0 <- c(1,1,1)
y <- rnorm(n,mean=mus0[cp0],sd=sqrt(sds0[cp0]))
hist(y,breaks = 30,col="blue1")
#B. initial value####################################################
#prior for theta##############################
ep <- 0.5
s0 <- 0.1
c <- 0.5#prior for v of the beta distribution##############################
Nc <- n#initial value for the total number of components###############
delta <- c(1:n)#indicator variable############
s <- exp(-delta)#deterministic decreasing sequence
theta <-  matrix(1, nrow=2, ncol=n)#initial for theta###
v <- NULL
w <- NULL
u <- NULL#latent variable 1###################
Ncy <- NULL#number of components of each y
Nc_all <- NULL#record Nc
delta_all <- NULL#record delta
m_all <- NULL #record non-empty groups
#big loop##################################################################
I <- 1#total number of iterations
#big loop##################################################################
IN <- 5000
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
  ksi <- NULL
  d <- NULL
  nb <- 10000#number of iterations in step b
  mu <- matrix(nrow=nb,ncol=Nc)#matrix with column for components and row for iteration.
  lambda <- mu###final estimation mu,lambda
  mu[1, ] <- theta[1,1:Nc]
  lambda[1,] <- theta[2,1:Nc]# first row for initial value
  for (j in 1:Nc)
  {
    if (m[j]!=0){
      ksi[j] <- sum(y[oc[[j]]])  
      lambda1 <- ep+m[j]/2
      for (ib in 2:nb)
      {
        mu1 <- ksi[j]*lambda[ib-1,j]/(m[j]*lambda[ib-1,j]+s0)
        mu2 <- 1/(m[j]*lambda[ib-1,j]+s0)#sd
        mu[ib,j] <- rnorm(1,mean=mu1,sd=mu2)
        d <- crossprod(y[oc[[j]]]-mu[ib,j])
        lambda2 <- ep+d/2
        lambda[ib,j] <- rgamma(1,shape=lambda1,rate=lambda2)#1/sd
      }
    } else {
      mu[nb,j] <- rnorm(1,mean=0,sd=1/s0)
      lambda[nb,j] <- rgamma(1,shape=ep,rate=ep)
    }
  }
  theta[,1:Nc] <- rbind(mu[nb,],lambda[nb,])#results from last iteration
  
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
  p<-matrix(nrow=n,ncol=Nc)
  for (ie in 1:n){
    for (je in (1:Nc)){
      if (s[je]>u[ie]){
        p[ie,je]<-w[je]*dnorm(y[ie],mean=theta[1,je],sd=sqrt(1/theta[2,je]))/s[je]
      }else{
        p[ie,je]<-0
      }
    }
    delta[ie]<-which(p[ie,]==max(p[ie,]))
    #delta[ie]<-sample(c(1:Nc),size=1,replace=FALSE, prob=p[ie,]/(sum(p[ie,])))
  }
  delta_all <- rbind(delta_all, delta)#record all group allocation

  I<-I+1
  if (I>IN)
    break
}


m_all[4000:5000]

cluster<-m_all[-1:-1000]
average<-NULL
for (i in 1:length(cluster)){
  average[i]<-mean(cluster[1:i])
}
average

f<-4716
unique(delta_all[f,])
which(delta_all[f,]==1)
which(delta_all[f,]==2)
sum(delta_all[f,]-delta_all[f-1,])




rho_k <- as.numeric(unname(unlist(acf(m_all,plot=FALSE))))[-1]
c_k <- which(rho_k<2/sqrt(IN))[1]
tao_k <- 0.5+sum(rho_k[1:(c_k-1)])
