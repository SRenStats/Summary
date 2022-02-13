setwd("C:/Users/DELL/Desktop/Role analysis")
require(ergm)
K=1 #distance for defining neighbourhood in start.R
MINSIZE=3 #minmum size of network
START="READ"
G=3 # number of groups / egonetwork types
# read in the network and attribute data
net<-network(as.matrix(read.table("ELadv.dat")))
atts<-read.csv("ELattr.csv");for (i in 1:ncol(atts)) set.vertex.attribute(net,names(atts)[i],atts[,i])
source("start.R")
source("likelihoods.R")
ego.terms<-c("edges", "mutual", "gwidegree(decay=0.8,fixed=TRUE)", "gwodegree(decay=0.8,fixed=TRUE)")
ergm.offset<-rep(0,N)
for (i in 1:N)
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."
source("initialise.R")
cat(date(), "starting intitialisation\n")
init<-init.egoergm(ego.terms, G) 
cat(date(), "finished intitialisation\n")

STEPS=100
cat(date(), "starting mixture model\n")
source("terms_ego_ergm.R")
source("mixtures_of_ergms.R")
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
cat(date(), "finished mixture model\n")
lambda<-out$lambda # membership results
group.theta<-out$theta # parameter results
EE.BIC<-out$EE.BIC # BIC result

require(blockmodeling)
cat(date(), "starting Regular Equivalence\n")
D<-REGE.for(M=as.sociomatrix(net))$E
clu<-cutree(hclust(d=as.dist(1-D),method="ward"),k=G) # other methods give similar results
cat(date(), "finished Regular Equivalence\n")
cat(date(), "starting 2-colouring\n")
source("maxmincoloring.R")
cat(date(), "finished 2-colouring\n")
cat(date(), "starting 3-colouring\n")
source("coloring3.R") 
cat(date(), "finished 3-colouring\n")
detach("package:igraph")
require(rapport)
library(rapport)
library(rapportools)

allstats<-cbind(get.vertex.attribute(net,"status"),clu,R2,R3,apply(init$lambda,1,which.max),apply(lambda,1,which.max))
lambdat<-matrix(1,6,6)
for (i in 1:6)
  for (j in 1:6) 
    lambdat[i,j]<-lambda.test(table(allstats[,i],allstats[,j]),direction=2)
colnames(lambdat)<-rownames(lambdat)<-c("reported", "Reg. Equiv.", "2-colour", "3-colour", "2-stage", "mixture")
print(lambdat)

# The match is for status==1 to be indicated by cluster 2 and status==2 by cluster 1
# Soft clustering looks useful because the mis-classified actors always have more uncertainty. 
mean(lambda[(get.vertex.attribute(net, "status")==1 & apply(lambda,1,which.max)==1),1]) # 6 actors
mean(lambda[(get.vertex.attribute(net, "status")==2 & apply(lambda,1,which.max)==1),1]) # 20 actors
# smaller cluster has lower mean lambda value

mean(lambda[(get.vertex.attribute(net, "status")==1 & apply(lambda,1,which.max)==2),2]) # 28 actors
mean(lambda[(get.vertex.attribute(net, "status")==2 & apply(lambda,1,which.max)==2),2]) # 7 actors
# smaller cluster has lower mean lambda value
save.image("lawyers_allmodels.RData")
