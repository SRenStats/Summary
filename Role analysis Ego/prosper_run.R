K=2 # how many hops out from the ego to include
MINSIZE=3 # minimum size of an egonetwork to consider
START="READ"
G=5 # number of egonetwork types; lenders and borrowers only
load("prosper_2010.rdata")
require(network)

source("start.R")
source("likelihoods.R")
ego.terms<-c("edges", "mutual", "gwidegree(decay=0.8,fixed=TRUE)", "gwodegree(decay=0.8,fixed=TRUE)")
ergm.offset<-rep(0,N)
for (i in 1:N)
  ergm.offset[i]<- -log(network.size(x[[i]]))*network.edgecount(x[[i]])
source("initialise.R")
cat(date(), "starting intitialisation\n")
init<-init.egoergm(ego.terms, G) 
cat(date(), "finished intitialisation\n")

source("terms_ego_ergm.R") # calculate the network statistics and create the formula

STEPS=50 # maximum EM iterations to do
cat(date(), "starting mixture model\n")
source("mixtures_of_ergms.R") # function to perform EM
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
cat(date(), "finished mixture model\n")
lambda<-out$lambda # store the mixture membership results of the EM algorithm
group.theta<-out$theta # store the parameter results of the EM algorithm
EE.BIC<-out$EE.BIC # store the BIC results of the EM algorithm

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
allstats<-cbind(get.vertex.attribute(net,"roles"),clu,R2,R3,apply(init$lambda,1,which.max),apply(lambda,1,which.max))
lambdat<-matrix(1,6,6)
for (i in 1:6)
  for (j in 1:6) 
    lambdat[i,j]<-lambda.test(table(allstats[,i],allstats[,j]),direction=2) # j predicts i
colnames(lambdat)<-rownames(lambdat)<-c("reported", "Reg. Equiv.", "2-colour", "3-colour", "2-stage", "mixture")
print(lambdat)

# look at soft clustering
roles<-unique(get.vertex.attribute(net,"roles"))[1:11]
# first col is correct maximally assigned actors and second are next highest but wrongly assigned
softclusmeans<-matrix(0,length(roles),2)
z<-apply(lambda,1,which.max)
alltab<-table(get.vertex.attribute(net,"roles"),z)
alltab<-t(t(alltab)/colSums(alltab))
for (r in 1:length(roles)) 
  {
  # first find the group that best matches
  tmp<-alltab[rownames(alltab)==roles[r],]
  cmax<-which.max(tmp)
  softclusmeans[r,1]=mean(lambda[(get.vertex.attribute(net,"roles")==roles[r]&z==cmax),cmax],na.rm=1)
  tmp[cmax]=0
  if (sum(tmp)>0)
    {
    cnext<-which.max(tmp)
    softclusmeans[r,2]=mean(lambda[(get.vertex.attribute(net,"roles")==roles[r]&z==cnext),cnext],na.rm=1)
    }
  }
# the maximally assigned actors have far higher mean lambda values showing that high uncertainty is
# strongly correlated with wrong maximal assignment.
save.image("prosper_allmodels.RData")
