n.mstepfunction2 <- function(group.theta.g,obs.S,label,g){
  a <- which (label==g)
  ans1 <- NULL
  for (i in 1:length(a)){
    ans1[i] <- pseudo.loglikelihood(group.theta.g,obs.S[[a[i]]])
  }
  ans<--sum(ans1)
  return(ans)
}

fit.mix.egoergm<-function(ego.terms, init, obs.S, G, p.ego.terms=NULL, p.group.theta=NULL)
  {
  #Nterms<-length(ego.terms)#dimension
  #ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
  #form<-nonsimp_update.formula(as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
  label<-initial.label
  TAU<-table(label)/N #colMeans(lambda)
  LL<-NaN
  LL_small <- NULL
  Nterms <- length(ego.terms)
  group.theta<-matrix(0,nrow=Nterms,ncol=G)
  
  lambda<-matrix(nrow=N,ncol=G)
  cat("iterations\tlog-likelood\n")
  ############# E-M iterations ##############
  for (l in 1:STEPS)    {
    oldLL<-LL
    #cat(l, " ", sep="")
    for (g in 1:G){
      fit <- nlm(n.mstepfunction2,group.theta[,g],obs.S,label,g)
      group.theta[,g]<-fit$estimate
      LL_small[g] <- -fit$minimum
    }
    LL<-sum(LL_small)

    # E-step
    #cat("E-step", l, "\t")
   # old_lambda<-lambda
    for (i in 1:N)    {
      for (g in 1:G) {
        lambda[i,g]<-log(TAU[g]) + pseudo.loglikelihood(group.theta[,g],obs.S[[i]]) }
      lambda[i,]<-lambda[i,]-max(lambda[i,])
    }
    lambda<-exp(lambda)
    lambda<-lambda/apply(lambda,1,sum)
    lambda[is.na(lambda)]<-1/G
    tmp<-apply(lambda,1,sum)
    lambda<-lambda/tmp # normalise lambda
    lambda[lambda==0]<-1e-8; lambda<-lambda/tmp
    #Update Label
    label<-apply(lambda,1,which.max)
    #Update weight
    TAU<-table(label)/N
    #cat("M-step", l, "\t")

    cat(l, "\t", sprintf("%.10f",LL),"\n")
    #if (max(abs(group.theta-old_group.theta))<tol)
    if (l>1)
      if ((LL-oldLL) < tol)
        {
        cat("Converged after ", l, "iterations\n")
        break
        }
    }
  #LL<-Mstepfunction(group.theta,obs.S,N,lambda,TAU,G)
  # higher is better
  EE.BIC<- 2*LL - (2*G*Nterms+G-1)*log(N) # usual BIC; higher is better
  return(list(EE.BIC=EE.BIC, theta=group.theta, label=label, G=G))
}


