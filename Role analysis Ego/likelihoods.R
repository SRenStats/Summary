######## Specify the pseudologlikelihood ########
pseudo.loglikelihood<-function(S,tmp.theta) 
  {
  loglike=sum(dbinom(S$response*S$weights,S$weights,
              inv.logit(S$offset+as.matrix(S$predictor)%*%tmp.theta),log=TRUE),na.rm=1)
  if(!is.finite(loglike)|loglike<LOWESTLL) 
    loglike=LOWESTLL# avoids numerical errors
  return(loglike)
  }
# returns the negative so that optimization is towards maximum
n.pseudo.loglikelihood<-function(S,tmp.theta) 
  -pseudo.loglikelihood(S,tmp.theta)

approx.loglikelihood<-function(S,tmp.theta,old.theta,M,form,ll0) pseudo.loglikelihood(S,tmp.theta)
  
# loglikelihood summed across all groups
Mstepfunction<-function(tmp.theta,S,N,lambda,TAU,G)
  {
  tmp.theta<-matrix(tmp.theta,nrow=G)
  ans=0
  for (g in 1:nrow(tmp.theta))
    ans = ans + mstepfunction(tmp.theta[g,],S,N,lambda[,g],TAU[g])
  return(ans)
  }
n.approx.loglikelihood<-function(S,tmp.theta,old.theta,M,form,ll0)
  -approx.loglikelihood(S,tmp.theta,old.theta,M,form,ll0)
mstepfunction<-function(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0)
  sum(lambda * (log(TAU) + sapply(S,approx.loglikelihood,tmp.theta,old.theta,M,form,ll0)))
n.mstepfunction<-function(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0)
  -mstepfunction(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0)

