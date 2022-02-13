init.egoergm<-function(ego.terms, G, p.ego.terms=NULL)
  {
  Nterms<-length(ego.terms)
  ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
  form<-ergm.update.formula(as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)

  ######### fit all ego-networks #############
  if (is.null(p.ego.terms))
    {
    theta<-matrix(0,N,Nterms)
    for (i in 1:N)
      theta[i,]<-ergmMPLE(form,output="fit")$coef
    # next couple of lines very ad-hoc but not an issue post EM convergence.
    theta[is.na(theta)]<-0
    theta[theta==-Inf]<- -1e6
    theta[theta==Inf]<- 1e6
    ############# initialisation ##############
    # k-means Clustering start
    if (G>1)
      {
      initial.clusters<-kmeans(theta, G, centers=apply(theta, 2, tapply,cutree(hclust(dist(theta)),G), mean),nstart=5)
      group.theta<-initial.clusters$centers
      group.theta<-matrix(group.theta,G,Nterms)
      }
    if (G==1)
      group.theta<-matrix(apply(theta,2,mean),nrow=1)
    lambda<-matrix(NaN,N,G)
    for (i in 1:N)
      for (g in 1:G)
        lambda[i,g]<-1/(sqrt(t(group.theta[g,]-theta[i,])%*%(group.theta[g,]-theta[i,]))+tol) # nugget to avoid zero
    lambda<-lambda/apply(lambda,1,sum) # normalise lambda
    print("Finished kmeans initialisation")
    }
  return(list(theta=theta, group.theta=group.theta, lambda=lambda, G=G))
  }
