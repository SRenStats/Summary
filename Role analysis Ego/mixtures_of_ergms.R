fit.mix.egoergm<-function(ego.terms, init, obs.S, G, p.ego.terms=NULL, p.group.theta=NULL)
  {
  Nterms<-length(ego.terms)
  ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
  form<-ergm.update.formula(as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
  lambda<-init$lambda
  group.theta<-init$group.theta
  TAU<-apply(lambda,2,sum)/N
  LL<-NaN
  cat("iterations\tlog-likelood\n")
  ############# E-M iterations ##############
  for (l in 1:STEPS)
    {
    oldLL<-LL
    #cat(l, " ", sep="")
    # E-step
    #cat("E-step", l, "\t")
    old_lambda<-lambda
    for (i in 1:N)
      {
      for (g in 1:G)
        lambda[i,g]<-log(TAU[g]) + approx.loglikelihood(obs.S[[i]],group.theta[g,],group.theta[g,]*0,M,form,0) 
      lambda[i,]<-lambda[i,]-max(lambda[i,])
      }
    lambda<-exp(lambda)
    lambda<-lambda/apply(lambda,1,sum)
    lambda[is.na(lambda)]<-1/G
    tmp<-apply(lambda,1,sum)
    lambda<-lambda/tmp # normalise lambda
    lambda[lambda==0]<-1e-8; lambda<-lambda/tmp
    TAU<-apply(lambda,2,sum)/N
    # M-step
    #cat("M-step", l, "\t")
    LL<-Mstepfunction(group.theta,obs.S,N,lambda,TAU,G)
    cat(l, "\t", sprintf("%.10f",LL),"\n")
    old_group.theta<-group.theta 
    for (g in 1:G)
      group.theta[g,]<-nlm(n.mstepfunction,group.theta[g,],S=obs.S,N=N,lambda=lambda[,g],TAU=TAU[g],group.theta[g,]*0,M,form,0)$estimate
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
  return(list(EE.BIC=EE.BIC, theta=group.theta, lambda=lambda, G=G))
  }
