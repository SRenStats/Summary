setwd("Z:/Role analysis Trade")
K=1#distance
MINSIZE=1#acceptable network size
START="SIM" 
#no definition of M,G  when repeating simulaitons
#G=3
#M=c(10,40,50) # M egonets per group, different weight
M<-c(50,50,50) #origional, same weight
N<-sum(M)  
source("start.R")#simulate egonets



ego.terms<-sim.ego.terms
ergm.offset<-rep(0,N)
for (i in 1:N)  
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."
source("initialise.R")#k-means initialization
#source("initialise20200113.R")#k-means initialization
init<-init.egoergm(ego.terms, G)
label<-init$label
#source("likelihoods20200113.R")#define functions
source("likelihoods.R")#define functions

source("terms_ego_ergm.R")#calculate network ststistics


initial.label<-sample(1:G,N,replace=TRUE)
STEPS=50
#source("mixtures_of_ergms20200113.R")
source("mixtures_of_ergms.R")
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm

table(sim.K,initial.label)
table(sim.K,out$label)
table(sim.K,apply(out$lambda,1,which.max))


