########1. pseudo.loglikelihood######
##(Specify the pseudologlikelihood) 
pseudo.loglikelihood<-function(tmp.theta,S) 
{
  loglike=sum(dbinom(S$response*S$weights,S$weights,
                     inv.logit(S$offset+as.matrix(S$predictor)%*%tmp.theta),log=TRUE),na.rm=1)#same as not using weight, either 0/full happens
  if(!is.finite(loglike)|loglike<LOWESTLL) #| or
    loglike=LOWESTLL# avoids numerical errors
  return(loglike)
}
## returns the negative so that optimization is towards maximum
n.pseudo.loglikelihood<-function(tmp.theta, S) {
  -pseudo.loglikelihood(S,tmp.theta)}

