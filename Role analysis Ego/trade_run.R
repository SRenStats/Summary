setwd("Z:/Role analysis Trade")
library(network)
K=1
MINSIZE=3
START="READ"
G=6 # number of groups / egonetwork types
# read all the data of the network and attribute data
#2014
year<-14
trade2<-as.matrix(read.csv("trade61.csv",header=TRUE))[(year*60-59):(year*60),]
gdp_value<-as.numeric(read.csv("gdp.csv",header=TRUE)[,(year+2)])/10000
#2016
# trade2<-as.matrix(read.csv("trade61.csv",header=TRUE))[(16*60-59):(16*60),]
# gdp_value<-as.numeric(read.csv("gdp.csv",header=TRUE)[,18])/10000
#dis<-as.matrix(read.csv("dis.csv",header=FALSE))
Covariate<-"GDP"
#Covariate<-NULL
trade<-ifelse(trade2>3*10^6,1,0)#binary matrix
#ergm(net~edges+gwidegree(decay=0.8,fixed=TRUE)+gwodegree(decay=0.8,fixed=TRUE))
#ego.terms<-c('edges', 'mutual', 'nodeicov("gdp")','nodeocov("gdp")')
ego.terms<-c('edges', 'mutual')
net<-as.network(as.matrix(trade))
network.density(net)
set.vertex.attribute(net,"gdp",gdp_value)
#set.edge.value(net,"distance",dis)#this works like weight rather than dyad
#set.edge.value(net,"amount",trade2016)


source("start.R")
network.density(net)
ergm.offset<-rep(0,N)
for (i in 1:N)
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."
#exp(-ergm.offset)
#i<-48
#plot(x[[i]])
den<-matrix(0,N,4)
for (i in 1:N){
  den[i,]<-c(network.size(x[[i]]),network.density(x[[i]]),network.dyadcount(x[[i]]),network.edgecount(x[[i]]))
}
mean(den[,2])




source("initialise.R")
cat(date(), "starting intitialisation\n")
init<-init.egoergm(ego.terms, G) 
cat(date(), "finished intitialisation\n")
table(apply(init$lambda,1,which.max))
init$group.theta
# network.vertex.names(net)[apply(init$lambda,1,which.max)==1]
# network.vertex.names(net)[apply(init$lambda,1,which.max)==2]
# network.vertex.names(net)[apply(init$lambda,1,which.max)==3]
# network.vertex.names(net)[apply(init$lambda,1,which.max)==4]


STEPS=100
source("likelihoods.R")
source("terms_ego_ergm.R") #calculate obs.S (MPLE ststistics) Y,X,W,N
source("mixtures_of_ergms.R")
cat(date(), "starting mixture model\n")
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
cat(date(), "finished mixture model\n")
lambda<-out$lambda # membership results
group.theta<-out$theta # parameter results
EE.BIC<-out$EE.BIC # BIC result
# The match is for status==1 to be indicated by cluster 2 and status==2 by cluster 1
# Soft clustering looks useful because the mis-classified actors always have more uncertainty. 
EE.BIC
table(apply(lambda,1,which.max))
group.theta


network.vertex.names(net)[apply(lambda,1,which.max)==1]
network.vertex.names(net)[apply(lambda,1,which.max)==2]
network.vertex.names(net)[apply(lambda,1,which.max)==3]
network.vertex.names(net)[apply(lambda,1,which.max)==4]


table(apply(init$lambda,1,which.max),apply(lambda,1,which.max))

mean(lambda[apply(lambda,1,which.max)==1,1]) # 6 actors













##########Compare with other methods###############
require(blockmodeling)
cat(date(), "starting Regular Equivalence\n")
D<-REGE.for(M=as.sociomatrix(net))$E
clu<-cutree(hclust(d=as.dist(1-D),method="ward.D"),k=G) # other methods give similar results
cat(date(), "finished Regular Equivalence\n")

cat(date(), "starting 2-colouring\n")
source("maxmincoloring.R")
cat(date(), "finished 2-colouring\n")

cat(date(), "starting 3-colouring\n")
source("coloring3.R") 
cat(date(), "finished 3-colouring\n")
detach("package:igraph")

library(rapport)
library(rapportools)
allstats<-cbind(clu,R2,R3,apply(init$lambda,1,which.max),apply(lambda,1,which.max))
lambdat<-matrix(1,5,5)
for (i in 1:5)
  for (j in 1:5) 
    lambdat[i,j]<-lambda.test(table(allstats[,i],allstats[,j]),direction=2)
colnames(lambdat)<-rownames(lambdat)<-c("Reg. Equiv.", "2-colour", "3-colour", "2-stage", "mixture")
print(lambdat)#the bigger the better


