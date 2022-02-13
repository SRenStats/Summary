setwd("Z:/Role analysis Trade")
library(network)
K=1
MINSIZE=2
START="READ"
G=3 # number of groups / egonetwork types
# read in the network and attribute data
net<-as.network(as.matrix(read.table("ELadv.dat")))
atts<-read.csv("ELattr.csv")
Covariate<-NULL

for (i in 1:ncol(atts)) {set.vertex.attribute(net,names(atts)[i],atts[,i])}



source("start2.R")


plot(x[[3]])

den_l<-matrix(0,N,4)
for (i in 1:N){
  den_l[i,]<-c(network.size(x[[i]]),network.density(x[[i]]),network.dyadcount(x[[i]]),network.edgecount(x[[i]]))
}



ergm.offset<-rep(0,N)
for (i in 1:N)
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."
ego.terms<-c("edges", "mutual", "gwidegree(decay=0.8,fixed=TRUE)", "gwodegree(decay=0.8,fixed=TRUE)")
 
#ego.terms<-c("edges", "mutual", "gwodegree(decay=0.8,fixed=TRUE)")


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
table(get.vertex.attribute(net, "status"),apply(lambda,1,which.max))
table(get.vertex.attribute(net, "status"),apply(init$lambda,1,which.max))
# The match is for status==1 to be indicated by cluster 2 and status==2 by cluster 1
# Soft clustering looks useful because the mis-classified actors always have more uncertainty. 
mean(lambda[(get.vertex.attribute(net, "status")==1 & apply(lambda,1,which.max)==1),1]) # 6 actors
mean(lambda[(get.vertex.attribute(net, "status")==2 & apply(lambda,1,which.max)==1),1]) # 20 actors
# smaller cluster has lower mean lambda value

mean(lambda[(get.vertex.attribute(net, "status")==1 & apply(lambda,1,which.max)==2),2]) # 28 actors
mean(lambda[(get.vertex.attribute(net, "status")==2 & apply(lambda,1,which.max)==2),2]) # 7 actors
# smaller cluster has lower mean lambda value
#save.image("lawyers_allmodels.RData")

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
allstats<-cbind(get.vertex.attribute(net,"status"),clu,R2,R3,apply(init$lambda,1,which.max),apply(lambda,1,which.max))
lambdat<-matrix(1,6,6)
for (i in 1:6)
  for (j in 1:6) 
    lambdat[i,j]<-lambda.test(table(allstats[,i],allstats[,j]),direction=2)
colnames(lambdat)<-rownames(lambdat)<-c("reported", "Reg. Equiv.", "2-colour", "3-colour", "2-stage", "mixture")
print(lambdat)#the bigger the better


