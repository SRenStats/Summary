setwd("~/Role analysis Trade")
require(boot) 
require(ergm) 
require(sna) 
LOWESTLL=-1e8 # lowest allowed loglikelihood to avoid numerical instabilities
tol=1e-6 # for convergence checking
# K=1#distance
# MINSIZE=1#acceptable network size
# START="SIM" 
#no definition of M,G  when repeating simulaitons
G=3
#M=c(10,40,50) # M egonets per group, different weight
M<-c(20,20,20) #origional, same weight
N<-sum(M)  #three groups with same number of nodes within each group
###################################################
###1. Simulate networks##########
require(statnet)
### create the underlying parameters
sim.theta<-rbind(c(-3,0,1), c(-1,-2,-1), c(-2,0,2)) # easy version
#sim.theta<-rbind(c(-0.05,0.05,0), c(-0.05,-0.05,-0.05), c(-0.05,0,0.05)) # hard version
sim.ego.terms<-c("edges", "gwesp(0.8,fixed=TRUE)", "gwdegree(decay=0.8,fixed=TRUE)")
#G=nrow(sim.theta)
### assign the egos to groups
sim.K<-NULL
for (g in 1:G)
  sim.K<-c(sim.K,rep(g,M[g]))
### simulate the egos
#NN=20 # number of nodes per group
ergmformula <- paste("~", paste(sim.ego.terms,collapse="+"),sep="")
sim.x<-list()
for (i in 1:G)
  sim.x[[i]]<-simulate(as.formula(paste("network(M[i],directed=FALSE)",ergmformula)), 
                       coef=c(sim.theta[(1:G)[sim.K[i]],]))

#save(sim.ego.terms,sim.theta,sim.K,sim.x, file="simulated_networks.rdata")
eps<-0.05
x11<-as.matrix.network(sim.x[[1]])
x22<-as.matrix.network(sim.x[[2]])
x33<-as.matrix.network(sim.x[[3]])
x21<-matrix(rbinom(M[2]*M[1],1,eps),nrow=M[2],ncol=M[1])
x31<-matrix(rbinom(M[3]*M[1],1,eps),nrow=M[3],ncol=M[1])
x32<-matrix(rbinom(M[3]*M[2],1,eps),nrow=M[3],ncol=M[2])
x1<-cbind(x11,t(x21),t(x31))
x2<-cbind(x21,x22,t(x32))
x3<-cbind(x31,x32,x33)
x<-rbind(x1,x2,x3)#x is a symmetric matrix

# library(igraph)
# graph<-graph_from_adjacency_matrix(x)
# V(graph)$color<-sim.K
# edge.start<-ends(graph,es=E(graph),names=F)
# edge.col <- V(graph)$color[edge.start]
# plot(graph, layout=layout_with_kk(graph),vertex.label=NA,vertex.size=4,edge.color=edge.col, edge.curved=.1,edge.arrow.size=0.5,edge.arrow.width=0.3)  
rm(sim.x,x1,x2,x3,x11,x22,x33,x21,x31,x32)
source("likelihoods3.R")
ego.terms<-sim.ego.terms
Nterms<-length(ego.terms)#dimension


#initialisation for label
max.iter<-10
label<-sample(1:G,N,replace = TRUE)


iter <- 0
RUN <- T
llh_old<-NA

while (RUN){
  #weight estimation and current likelihood
  TAU<-table(label)/N
  
  
  #theta estimation
  y<-list()
  group.theta<-matrix(0,G,Nterms)
  #ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
  form<-nonsimp_update.formula(as.formula(paste("y[[i]]",ergmformula)),y[[i]] ~ .)
  loglike<-NULL
  for (i in 1:G){
    y[[i]]<-network(x[which(label==i),which(label==i)],directed = FALSE)
    # ergm(form,estimate="MPLE") #originional version
    data1 <- ergmMPLE(form)#get the Y~X form
    model1<-glm(response~predictor[,-1],weights=weights,data=data1,family = "binomial")
    group.theta[i,]<-model1$coefficients
    loglike[i]<-logLik(model1) } #log-likelihood of this group
  llh_new<-sum(loglike)
  
  
  
  #update label E-STEP
  lambda<-matrix(nrow=N,ncol=G)
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
  
  label<-apply(lambda,1,which.max)
  
  if(!is.na(llh_old)){
    if ((llh_new-llh_old)/llh_old<tol)
      RUN <- F
  }
  
  
  llh_old<-llh_new
  
  iter= iter+1
  if(iter >= max.iter){
    print("Maximum Iterations Achieved: Fit Suspect")
    RUN = F #STOP interation
  }
}
table(sim.K,label)











