require(boot) 
require(ergm) 
require(sna) 
LOWESTLL=-1e8 # lowest allowed loglikelihood to avoid numerical instabilities
tol=1e-6 # for convergence checking


#fit.mix.egoergm<-function(ego.terms, init, obs.S, G, p.ego.terms=NULL, p.group.theta=NULL)
#  {

  #ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
  #form<-nonsimp_update.formula(as.formula(paste("y[[i]]",ergmformula)),y[[i]] ~ .)
  #lambda<-init$lambda
  #group.theta<-init$group.theta
  #label<-init$label

  # obs.S<-list()#list of every ego-network with 4 objects
  # obs.S_o<-list()
  # y<-list()
  # y_o<-list()
  LL<-NULL
  lambda<-matrix(nrow=N,ncol=G)
  #group.lambda<-rep(NA,G)
  cat("iterations\tlog-likelood\n")
  ############# E-M iterations ##############
  for (l in 1:STEPS)    {
    oldLL<-LL
    #cat(l, " ", sep="")
    # E-step
    #cat("E-step", l, "\t")
   # old_lambda<-lambda

    for (i in 1:N)      {
      for (g in 1:G)      {
      lambda[i,g]<-log(TAU[g]) + row.log.likelihood(x,group.theta[g,],i,g,label)    
        }
      lambda[i,]<-lambda[i,]-max(lambda[i,])
      }
    lambda<-exp(lambda)
    lambda<-lambda/apply(lambda,1,sum)
    lambda[is.na(lambda)]<-1/G
    tmp<-apply(lambda,1,sum)
    lambda<-lambda/tmp # normalise lambda
    lambda[lambda==0]<-1e-8; lambda<-lambda/tmp
    #Update weight
    TAU<-apply(lambda,2,sum)/N
    #The possible problem is with the label
    label<-apply(lambda,1,which.max)
    # M-step
    #cat("M-step", l, "\t")
    La<-matrix(nrow=N,ncol=G)
    for (i in 1:N)      {
      for (g in 1:G)      {
        La[i,g]<-lambda[i,g]*(log(TAU[g]) + row.log.likelihood(x,group.theta[g,],i,g,label))          }
    }
    LL<-sum(La)
    cat(l, "\t", sprintf("%.10f",LL),"\n")
    old_group.theta<-group.theta 
    for (g in 1:G)
      group.theta[g,]<-nlm(n.mstepfunction2,group.theta[g,],g,lambda[,g],TAU[g],label,x)$estimate
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
#  return(list(EE.BIC=EE.BIC, theta=group.theta, lambda=lambda, label=label,G=G))
#  }
