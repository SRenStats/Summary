#########0.test u######
u<-(0:50)/50
theta1<-1-sin(2*u*pi)
theta2<-cos(3*u*pi)-1
plot(u,theta1,ylim = c(-2,2))
points(u,theta2)
######1.generate cluster membership vactor z and theta######
n<-100
T<-51
z<-rbinom(n,1,0.5)+1
u<-(0:50)/50
theta<-matrix(nrow=2,ncol=T)
theta[1,]<-0.5*cos(3*u*pi)
theta[2,]<--0.5*cos(3*u*pi)-0.75
#plot(u,theta[1,],ylim = c(-1.25,0.5))
#points(u,theta[2,])
######2.generate uptriangle network######
data<-array(0,dim=c(n,n,T))
for (t in 1:T)
{
    for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      pr0<-1/(1+exp(theta[1,t]+theta[2,t]))
      data[i,j,t]<-rbinom(1,1,1-pr0)
    }
  }
}
######3.#####