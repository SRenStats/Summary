########1. pseudo.loglikelihood######
##(Specify the pseudologlikelihood) 
pseudo.loglikelihood<-function(S,tmp.theta) 
{
  loglike=sum(dbinom(S$response*S$weights,S$weights,
                     inv.logit(S$offset+as.matrix(S$predictor)%*%tmp.theta),log=TRUE),na.rm=1)#same as not using weight, either 0/full happens
  if(!is.finite(loglike)|loglike<LOWESTLL) #| or
    loglike=LOWESTLL# avoids numerical errors
  return(loglike)
}
## returns the negative so that optimization is towards maximum
n.pseudo.loglikelihood<-function(S,tmp.theta) {
  -pseudo.loglikelihood(S,tmp.theta)}

###2. approx.loglikelihood #####
approx.loglikelihood<-function(S,tmp.theta,old.theta,M,form,ll0){
  pseudo.loglikelihood(S,tmp.theta)}

n.approx.loglikelihood<-function(S,tmp.theta,old.theta,M,form,ll0)
{  -approx.loglikelihood(S,tmp.theta,old.theta,M,form,ll0) }

####3. Mstepfunction#####
#loglikelihood summed across all groups
mstepfunction<-function(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0){
  sum(lambda * (log(TAU) + sapply(S,approx.loglikelihood,tmp.theta,old.theta,M,form,ll0)))}

n.mstepfunction<-function(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0){
  -mstepfunction(tmp.theta,S,N,lambda,TAU,old.theta,M,form,ll0)}

Mstepfunction<-function(tmp.theta,S,N,lambda,TAU,G)
{
  tmp.theta<-matrix(tmp.theta,nrow=G)
  ans=0
  for (g in 1:nrow(tmp.theta))
    ans = ans + mstepfunction(tmp.theta[g,],S,N,lambda[,g],TAU[g])
  return(ans)
}

#single g 
row.log.likelihood<-function(fullmatrix,group.theta.g,row,rowlabel,fulllabel){
  S<-NULL#list of every ego-network with 4 objects
  nodes<-unique(c(row,which(fulllabel==rowlabel)))
  offset<--log(length(nodes))
  y<-fullmatrix[nodes,nodes]
  for (i in 2:length(nodes)){
    y1<-y;y1[i,1]<-1;y1[1,i]<-1;y1<-network(y1,directed = FALSE)
    y0<-y;y0[i,1]<-0;y0[1,i]<-0;y0<-network(y0,directed = FALSE)
    S<-rbind(S,c(y[i,1],summary(as.formula(paste("y1",ergmformula)))-summary(as.formula(paste("y0",ergmformula)))))
  }
  S<-as.matrix(S)
  loglike<-sum(dbinom(S[,1],1,
                      inv.logit(offset+as.matrix(S[,-1]%*%group.theta.g)),log=TRUE),na.rm=1)#same as not using weight, either 0/full happens
  if(!is.finite(loglike)|loglike<LOWESTLL) #| or
    loglike=LOWESTLL# avoids numerical errors
  return(loglike)
}

mstepfunction2<-function(group.theta.g,g,lambda.g,TAU.g,fulllabel,fullmatrix)
  {
  a<-rep(NA,nrow(fullmatrix))
  for (i in 1:nrow(fullmatrix)){
a[i]<-lambda.g[i] * (log(TAU.g)[i] + row.log.likelihood(fullmatrix,group.theta.g,i,g,fulllabel))
  }
  return(sum(a))
 }

n.mstepfunction2<-function(group.theta.g,g,lambda.g,TAU.g,fulllabel,fullmatrix)
-mstepfunction2(group.theta.g,g,lambda.g,TAU.g,fulllabel,fullmatrix)

  