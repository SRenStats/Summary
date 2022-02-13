library(mcsm)
#demo(Chapter.7)
#Section 7.2, Energy example 7.3
data(Energy)
x <- Energy$Girls
n <- length(x)
nsim <- 10^4
a <- 3
b <- 3
tau2 <- 10
theta0 <- 5
xbar <- mean(x)
sh1 <- (n/2)+a
sigma <- rep(0,nsim)                  #init arrays
theta <- sigma
sigma[1] <- 1/rgamma(1,shape=a,rate=b)      #init chains
B <- sigma[1]/(sigma[1]+n*tau2)
theta[1] <- rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
 for (i in 2:nsim){
      B <- sigma[i-1]/(sigma[i-1]+n*tau2)
      theta[i] <- rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
      ra1 <- (1/2)*(sum((x-theta[i])^2))+b
      sigma[i] <- 1/rgamma(1,shape=sh1,rate=ra1)
 }

par(mfrow=c(1,1),mar=c(4,4,2,1))

plot(density(theta[-1:-5000]),type='l')

hist(log(theta[theta>0]), nclass=140, col="grey", freq=FALSE, main="",ylab="", xlab=expression(theta), 
     xlim=as.vector(quantile(log(theta[theta>0]), prob=c(.005,.995))))

hist(log(sigma), nclass=150, col="sienna", freq=FALSE, ,main="",ylab="", xlab=expression(sigma^2), 
     xlim=as.vector(quantile(log(sigma),prob=c(.005,.995))))


data(Energy)
x <- Energy$Girls
n <- length(x)
nsim <- 10^4
a <- 3
b <- 3
tau2 <- 10
theta0 <- 5
xbar <- mean(x)
sh1 <- (n/2)+a
sigma <- rep(0,nsim)                  #init arrays
theta <- sigma
sigma[1] <- 1/rgamma(1,shape=a,rate=b)      #init chains
B <- sigma[1]/(sigma[1]+n*tau2)
theta[1] <- rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
for (i in 2:nsim){
  B <- sigma[i-1]/(sigma[i-1]+n*tau2)
  theta[i] <- rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
  theta[i] <- B*theta0+(1-B)*xbar
  ra1 <- (1/2)*(sum((x-theta[i])^2))+b
  #sigma[i] <- 1/rgamma(1,shape=sh1,rate=ra1)
  sigma[i] <- ra1/sh1
}