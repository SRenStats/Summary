#1.generate data from a mixture model#########################################
set.seed(111)
n <- 50
nc0<-3
w0<- c(1/3,1/3,1/3)
cp0 <- sample(1:nc0,prob=w0,size=n,replace=TRUE)#component
mus0 <- c(-4,0,8)
sds0 <- c(1,1,1)
y <- rnorm(n,mean=mus0[cp0],sd=sqrt(sds0[cp0]))
#plot(density(y))
hist(y,breaks = 30,col="blue1")
#hist(y,breaks = c(0.5*(-12:4),3:5,0.5*(12:20)),col="blue1")
#2. initial value####################################################
#prior for theta##############################B
ep<-0.5
s<-0.1
a<-0.1 #prior for c
b<-0.1
c<-rgamma(1,shape=a,rate=b)
#initial
#n <- 50#number of obsevations
Nc <- n#######################################D
delta <- c(1:n)#indicator variable############E
w <- exp(-1*delta)#weight of each component###D
theta <-  matrix(1, nrow=2, ncol=n)#initial for theta###B
v <- rbeta(n,1,c)#initial for beta############C
u <- NULL#latent variable 1###################A
#big loop#################################################################
#big loop##################################################################
I<-1
IN<-30
repeat{
  #A.update u#######################################################
  for (ia in 1:n){
    u[ia] <- runif(1,0,w[delta[ia]])
  }
  #observations within each component

  oc<-vector("list",Nc)
  for (j in 1:Nc)
  {
    oc[[j]]<-which(delta==j)
      }
  #B.update theta=(mu,lambda)###################################################
  m<-NULL
  ksi <- NULL
  d<-NULL
  nb<- 5000#number of iterations in step b
  mu<-matrix(nrow=nb,ncol=Nc)#matrix with column for components and row for iteration.
  lambda<-mu
  mu[1, ] <- theta[1,1:Nc]
  lambda[1,] <- theta[2,1:Nc]# first row for initial value
  for (j in 1:Nc)
  {
    m[j] <-length(oc[[j]])
    if (m[j]!=0){
      ksi[j] <- sum(y[oc[[j]]])  
      lambda1 <- ep+m[j]/2
      for (ib in 2:nb)
      {
        mu1 <- ksi[j]*lambda[ib-1,j]/(m[j]*lambda[ib-1,j]+s)
        mu2 <- 1/(m[j]*lambda[ib-1,j]+s)#sd
        mu[ib,j] <- rnorm(1,mean=mu1,sd=sqrt(mu2))
        d <- crossprod(y[oc[[j]]]-mu[ib,j])
        lambda2 <- ep+d/2
        lambda[ib,j] <- rgamma(1,shape=lambda1,rate=lambda2)#1/sd
      }
    } else {
      mu[nb,j]<-rnorm(1,mean=0,sd=sqrt(1/s))
      lambda[nb,j]<-rgamma(1,shape=ep,rate=ep)
    }
  }
  theta<-rbind(mu[nb,],lambda[nb,])#results from last iteration
  #C.update v##########################################
  alpha<-NULL
  beta<-NULL
  for (j in 1:Nc)
  {
    alpha0 <- u[oc[[j]]]/prod(1-v[1:(j-1)])
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
    }#else vj keep unchaging
    #else{
    # v[j]<-rbeta(1,1,c)
    #}
  }  
  #D.update w and the number of components Nc###########################################
  u_star<-min(u)
  sum_w<-NULL
  w[1]<-v[1]
  sum_w[1]<-w[1]
  for (id in (2:Nc)){
    w[id]<-v[id]*prod(1-v[1:(id-1)])
    sum_w[id]<-sum(w[1:id])
  }
  Nc<-max(which(sum_w<1-u_star),delta)
  w<-w[1:Nc]
  #w<-w[order(w,decreasing = TRUE)[1:Nc]]
  v<-v[1:Nc]
  #E.update delta######################################################################
  p<-matrix(nrow=n,ncol=Nc)
  for (ie in 1:n){
    for (je in (1:Nc)){
      if (w[je]>u[ie]){
        p[ie,je]<-dnorm(y[ie],mean=theta[1,je],sd=sqrt(1/theta[2,je]))
      }else{
        p[ie,je]<-0
      }
    }
    delta[ie]<-which(p[ie,]==max(p[ie,]))
    #delta[ie]<-sample(c(1:Nc),size=1,replace=FALSE, prob=p[ie,]/(sum(p[ie,])))
  }
  #F.update c##########################################################################gibbs sampling
  c_g<-NULL
  eta<-NULL
  eta[1]<-rbeta(1,c+1,n)
  nf<-1000
  for (iff in 2:nf){
    pi<-(a+Nc-1)/(n*b-n*log(eta[iff-1])+a+Nc-1)
    u<-runif(1)
    if (u<pi){
      c_g[iff]<-rgamma(1,shape=a+Nc,rate=b-log(eta[iff-1]))
    }else{
      c_g[iff]<-rgamma(1,shape=a+Nc-1,rate=b-log(eta[iff-1]))
    }
    eta[iff]<-rbeta(1,c_g[iff]+1,n)
  }
  c<-c_g[nf]  

  I<-I+1
  if (I>IN)
    break
}





