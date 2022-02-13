################Get familiar with the data########3
library(igraph)
library(igraphdata)
library(intergraph)
#data(package="igraphdata")
data(karate)
plot(karate)

karate<-asNetwork(karate)
library(network)
library(ergm)
#Faction is the same as color
get.vertex.attribute(karate,"Faction")#actual group
get.vertex.attribute(karate,"label")[which(get.vertex.attribute(karate,"Faction")==1)]
get.edge.attribute(karate,"weight")# the number of common activities
summary(karate~edges+gwdegree(decay=0.8,fixed=TRUE))
#############Start #####################
setwd("Z:/Role analysis Trade")
K=1
MINSIZE=6
START="READ"
G=2 # number of groups / egonetwork types
# read in the network and attribute data
net<-karate
Covariate<-NULL
source("start.R")

den<-matrix(0,N,4)
for (i in 1:N){
  den[i,]<-c(network.size(x[[i]]),network.density(x[[i]]),network.dyadcount(x[[i]]),network.edgecount(x[[i]]))
}

ergm.offset<-rep(0,N)
for (i in 1:N)
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."
ego.terms<-c('edges','gwdegree(decay=0.8,fixed=TRUE)')
#ego.terms<-c('edges')


source("initialise.R")
cat(date(), "starting intitialisation\n")
init<-init.egoergm(ego.terms, G) 
cat(date(), "finished intitialisation\n")



STEPS=100
cat(date(), "starting mixture model\n")
source("likelihoods.R")
source("terms_ego_ergm.R") #calculate obs.S (MPLE ststistics) Y,X,W,N
source("mixtures_of_ergms.R")
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
cat(date(), "finished mixture model\n")
lambda<-out$lambda # membership results
group.theta<-out$theta # parameter results
EE.BIC<-out$EE.BIC # BIC result

table(apply(init$lambda,1,which.max))

table(apply(lambda,1,which.max))

table(get.vertex.attribute(net,"Faction"),apply(init$lambda,1,which.max))

table(get.vertex.attribute(net,"Faction"),apply(lambda,1,which.max))

get.vertex.attribute(net,"label")[apply(lambda,1,which.max)==2]#not even cluster two hubs


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
allstats<-cbind(get.vertex.attribute(net,"Faction"),clu,R2,R3,apply(init$lambda,1,which.max),apply(lambda,1,which.max))
lambdat<-matrix(1,6,6)
for (i in 1:6)
  for (j in 1:6) 
    lambdat[i,j]<-lambda.test(table(allstats[,i],allstats[,j]),direction=2)
colnames(lambdat)<-rownames(lambdat)<-c("reported", "Reg. Equiv.", "2-colour", "3-colour", "2-stage", "mixture")
print(lambdat)#the bigger the better


