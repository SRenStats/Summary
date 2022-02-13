require(boot) 
require(ergm) 
require(sna) 
LOWESTLL=-1e8 # lowest allowed loglikelihood to avoid numerical instabilities
tol=1e-6 # for convergence checking

init.egoergm<-function(ego.terms, G, p.ego.terms=NULL)
{
  Nterms<-length(ego.terms)
  #ego.terms<-c("edges", "mutual", "nodeicov("gdp")","nodeocov("gdp")")
  ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep='')
  #form<-ergm.update.formula(as.formula(paste("x[[i]]",ergmformula)),x[[i]] ~ .)
  form<-nonsimp_update.formula(as.formula(paste("y[[i]]",ergmformula)),y[[i]] ~ .)
  label<-sample(1:G,N,replace = TRUE)
  y<-list()
  ######### fit all ego-networks #############
  if (is.null(p.ego.terms))
  {
    group.theta<-matrix(0,G,Nterms)
    for (i in 1:G)
    {
      y[[i]]<-network(x[which(label==i),which(label==i)],directed = FALSE)
      group.theta[i,]<-ergmMPLE(form,output="fit")$coef    
      #theta[i,]<-ergm(form,estimate="MPLE")$coef
    }
    # next couple of lines very ad-hoc but not an issue post EM convergence.
    # theta[is.na(theta)]<-0
    # theta[theta==-Inf]<- -1e6
    # theta[theta==Inf]<- 1e6
    ############# initialisation ##############
    # theta<-matrix(0,N,Nterms)
    # ergm.offset<-rep(0,N)
    # for (i in 1:N){
    #   theta[i,]<-group.theta[label[i],]
    #   ergm.offset[i]<- size[label[i]]
    # }
    # 
    #lambda<-label
    print("Finished kmeans initialisation")
  }
  return(list(group.theta=group.theta, label=label))
}
