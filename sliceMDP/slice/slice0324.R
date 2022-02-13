#1.generate data from a mixture model#########################################
n <- 50
nc0<-3
w0<- c(1/3,1/3,1/3)
cp0 <- sample(1:nc0,prob=w0,size=n,replace=TRUE)#component
mus0 <- c(-8,0,8)
sds0 <- sqrt(c(1,1,1))
y <- rnorm(n,mean=mus0[cp0],sd=sds0[cp0])
plot(density(y))
#2. initial value####################################################
#prior for theta##############################B
ep<-0.5
s<-0.1
c<-1
#initial settings
Nc <- n#######################################E
delta <- c(1:n)#indicator variable############D
w <- exp(-1*delta)#weight of each component###E
#theta <-  matrix(1, nrow=2, ncol=n)#initial for theta###B
v <- rbeta(n,1,c)#initial for beta############C
u <- NULL#latent variable 1###################A
theta<-vector("list")
theta[[1]]<-matrix(1,nrow=2,ncol=Nc)
#big loop#################################################################
#big loop#################################################################
IN<-10000
for (I in (2:IN)){
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
  m <- NULL
  ksi <- NULL
  d <- NULL
  mu <- NULL
  lambda <- NULL# first row for initial value
  #theta<-vector("list")
  for (j in 1:Nc)
  {
    m[j] <-length(oc[[j]])
    if (m[j]!=0){
      ksi[j] <- sum(y[oc[[j]]])  
      lambda1 <- ep+m[j]/2
      mu[j]<-ksi[j]*theta[[I-1]][2,j]/(m[j]*theta[[I-1]][2,j]+s)
      d <- crossprod(y[oc[[j]]]-mu[j])
      lambda2 <- ep+d/2
      lambda[j] <- lambda1/lambda2#1/sd
    } else {
      mu[j] <- 0
      lambda[j]<-1
    }
  }
  theta[[I]]<-rbind(mu,lambda)
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
    else{
      v[j]<-rbeta(1,1,c)
    }
  }  
  #D.update delta######################################################################
  p<-matrix(nrow=n,ncol=Nc)
  for (id in 1:n){
    for (jd in (1:Nc)){
      if (w[jd]>u[id]){
        p[id,jd]<-dnorm(y[id],mean=theta[[I]][1,jd],sd=1/theta[[I]][2,jd])
      }else{
        p[id,jd]<-0
      }
    }
    delta[id]<-which(p[id,]==max(p[id,]))[1]
    #delta[id]<-sample(c(1:Nc),size=1,replace=FALSE, prob=p[id,]/(sum(p[id,])))
  }
  #E.update w and the number of components Nc###########################################
  u_star<-min(u)
  sum_w<-NULL
  w[1]<-v[1]
  sum_w[1]<-w[1]
  for (ie in (2:Nc)){
    w[ie]<-v[ie]*prod(1-v[1:(ie-1)])
    sum_w[ie]<-sum(w[1:ie])
  }
  Nc<-max(which(sum_w<1-u_star),delta)
  w<-w[order(w,decreasing = TRUE)[1:Nc]]
  #w<-w[1:Nc]
  v<-v[1:Nc]

}
