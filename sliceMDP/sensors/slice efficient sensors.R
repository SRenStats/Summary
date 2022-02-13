#Dirichlet process prior
#slice_efficient(beta_all,100,1)#1.data/2.number of iterations/3.c
#slice_efficient <- function(y,IN,c)

#y <- beta0#OLS estimates
y <- beta_training
c <- 3
library(mvtnorm)
n <- nrow(y)
p <- ncol(y)
#prior for theta##############################
mu0 <- rep(0,p)
k0 <- 1
v0 <- p
T0 <- diag(p)
#initial values###############
Nc <- n#initial value for number of clusters
delta <- c(1:n)#indicator variable for allocation############E
s <- exp(-delta)#deterministic decreasing sequence
theta <-  array(1, dim=c((p+1),p,n))#initial for theta###B
Ncy <- NULL
Nc_all <- NULL
delta_all <- NULL
m_all <- NULL
v <- NULL
w <- NULL
u <- NULL#latent variable 1
#big loop##################################################################
I <- 1#total number of iterations
IN <- 10000
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
      #lambda <- lambda1 * solve(lambda2)
      mu <- rmvnorm(1,mean = mu1, sigma = solve(mu2*lambda))#vector
      #mu <- mu1
    } else {
      lambda <- rWishart(1,df=v0,Sigma=T0)[,,1]
      #lambda <- v0 * solve(T0)
      mu <- rmvnorm(1,mean = mu0, sigma = solve(k0*lambda))
      #mu <- mu0
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
m_all[-1:-6000]
delta
length(which(delta==1))
length(which(delta==3))
colMeans(y[which(delta==1),])
colMeans(y[which(delta==3),])
w[unique(delta)]
theta[1,,unique(delta)]
