# code to re-run many simulated data examples to look at errors and biases of parameter estimates
R<-50
G<-3
alltheta<-array(NaN,c(R,G,3))
alllambda<-array(NaN,c(R,N,G))
allerrorcount<-rep(NaN,R)
initalltheta<-array(NaN,c(R,G,3))
initalllambda<-array(NaN,c(R,N,G))
initallerrorcount<-rep(NaN,R)
for (m in 1:R)
  {
  source("sim_run.R")
  #now address label switching problem
  perm<-apply(table(sim.K,apply(init$lambda,1,which.max)),1,which.max)
  if (length(perm==3))
    {
    init$lambda<-init$lambda[,perm]
    initalllambda[m,,]<-init$lambda
    init$group.theta<-init$group.theta[perm,]
    initalltheta[m,,]<-init$group.theta
    initallerrorcount[m]<-sum(apply(init$lambda,1,which.max)!=sim.K)
    }
  #now address label switching problem
  perm<-apply(table(sim.K,apply(lambda,1,which.max)),1,which.max)
  if (length(perm==3))
    {
    lambda<-lambda[,perm]
    alllambda[m,,]<-lambda
    group.theta<-group.theta[perm,]
    alltheta[m,,]<-group.theta
    allerrorcount[m]<-sum(apply(lambda,1,which.max)!=sim.K)
    }
  cat(m, "repetitions\n")
  }

