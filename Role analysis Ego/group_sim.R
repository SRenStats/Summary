setwd("Z:/Role analysis Trade")
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

library(igraph)
graph<-graph_from_adjacency_matrix(x)
V(graph)$color<-sim.K
edge.start<-ends(graph,es=E(graph),names=F)
edge.col <- V(graph)$color[edge.start]
plot(graph, layout=layout_with_kk(graph),vertex.label=NA,vertex.size=4,edge.color=edge.col, edge.curved=.1,edge.arrow.size=0.5,edge.arrow.width=0.3)  


rm(sim.x,x1,x2,x3,x11,x22,x33,x21,x31,x32)


ego.terms<-sim.ego.terms
Nterms<-length(ego.terms)#dimension
TAU<-c(0.33,0.33,0.34) #colMeans(lambda)
label<-sample(1:G,N,replace = TRUE)
y<-list()
group.theta<-matrix(0,G,Nterms)


#ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep="")
form<-nonsimp_update.formula(as.formula(paste("y[[i]]",ergmformula)),y[[i]] ~ .)

for (i in 1:G)
{
  y[[i]]<-network(x[which(label==i),which(label==i)],directed = FALSE)
  group.theta[i,]<-ergmMPLE(form,output="fit")$coef    
}





source("likelihoods3.R")#define functions
STEPS=5
source("mixtures_of_ergms3.R")



# out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
# lambda<-out$lambda
# group.theta<-out$theta
# EE.BIC<-out$EE.BIC


# table(sim.K,apply(init$lambda,1,which.max))
# table(sim.K,apply(lambda,1,which.max))
# chisq.test(table(sim.K,apply(init$lambda,1,which.max)))
# chisq.test(table(sim.K,apply(lambda,1,which.max)))
# lambda.test(table(sim.K,apply(init$lambda,1,which.max)))
# lambda.test(table(sim.K,apply(lambda,1,which.max)))
# 
# mean(1==apply(lambda,1,which.max))
# group.theta
# sim.theta
# colMeans(lambda)

