setwd("Z:/Role analysis Trade")
#install.packages("SRFM.ERGMuserterms_3.10.0.tar.gz", repos = NULL, type = "source", INSTALL_opts = "--no-multiarch")
library(statnet)
#library(SRFM.ERGMuserterms)
library(clues)
#source("SRFM_ERGM.R")
require(boot) 
require(ergm) 
require(sna) 
require(statnet)


G=3

M<-c(50,50,50) #origional, same weight
N<-sum(M)  #three groups with same number of nodes within each group
sim.theta<-cbind(c(-3,0,1), c(-1,-2,-1), c(-2,0,2)) # easy version
sim.K<-NULL
for (g in 1:G)
  sim.K<-c(sim.K,rep(g,M[g]))


###1. SIMULATE NERWTOKS###
 ego.terms<-c("edges", "gwesp(0.8,fixed=TRUE)", "gwidegree(decay=0.8,fixed=TRUE)"); ergmformula <- paste("~", paste(ego.terms,collapse="+"),sep='')
sim.x<-list()
for (i in 1:G)
  sim.x[[i]]<-simulate(as.formula(paste("network(M[i],directed=TRUE)",ergmformula)),
                       coef=c(sim.theta[(1:G)[sim.K[i]],]))

#save(sim.ego.terms,sim.theta,sim.K,sim.x, file="simulated_networks.rdata")
eps<-0.05
x11<-as.matrix.network(sim.x[[1]])
x12<-matrix(rbinom(M[1]*M[2],1,eps),nrow=M[1],ncol=M[2])
x13<-matrix(rbinom(M[1]*M[3],1,eps),nrow=M[1],ncol=M[3])

x21<-matrix(rbinom(M[2]*M[1],1,eps),nrow=M[2],ncol=M[1])
x22<-as.matrix.network(sim.x[[2]])
x23<-matrix(rbinom(M[2]*M[3],1,eps),nrow=M[2],ncol=M[3])

x31<-matrix(rbinom(M[3]*M[1],1,eps),nrow=M[3],ncol=M[1])
x32<-matrix(rbinom(M[3]*M[2],1,eps),nrow=M[3],ncol=M[2])
x33<-as.matrix.network(sim.x[[3]])

x1<-cbind(x11,x12,x13)
x2<-cbind(x21,x22,x23)
x3<-cbind(x31,x32,x33)
x<-rbind(x1,x2,x3)#x is a symmetric matrix

rm(sim.x,x1,x2,x3,x11,x22,x33,x21,x31,x32)

network<-network(x,directed=T)

###2. GENERATE EGO-NETWORKS###
K=1
MINSIZE=3
START="READ"
Covariate<-NULL
net<-network

source("start20200113.R")

###3. ESTIMATION###
ergm.offset<-rep(0,N)
for (i in 1:N)
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."


# G<-3
# source("initialise20200113.R")
# cat(date(), "starting intitialisation\n")
# init<-init.egoergm(ego.terms, G) 
# cat(date(), "finished intitialisation\n")

#initial.label<-init$label



source("terms_ego_ergm.R") #calculate obs.S (MPLE ststistics) Y,X,W,N

initial.label <- sample(1:G,N,replace = TRUE)
STEPS=10
cat(date(), "starting mixture model\n")
source("likelihoods20200113.R")

source("mixtures_of_ergms20200113.R")
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
cat(date(), "finished mixture model\n")

table(sim.K,initial.label)
table(sim.K,out$label)

