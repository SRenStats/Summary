setwd("Z:/Role analysis Trade")
K=1#distance
MINSIZE=1#acceptable network size
START="SIM" 
#no definition of M,G  when repeating simulaitons
#G=3
#M=c(10,40,50) # M egonets per group, different weight
M<-c(50,50,50) #origional, same weight
N<-sum(M)  
source("start2.R")#simulate egonets

source("likelihoods.R")#define functions

ego.terms<-sim.ego.terms
ergm.offset<-rep(0,N)
for (i in 1:N)  
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."
source("initialise.R")#k-means initialization
init<-init.egoergm(ego.terms, G) 


source("terms_ego_ergm.R")#calculate network ststistics
STEPS=50
source("mixtures_of_ergms.R")
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
lambda<-out$lambda
group.theta<-out$theta
EE.BIC<-out$EE.BIC


table(sim.K,apply(init$lambda,1,which.max))
table(sim.K,apply(lambda,1,which.max))
# chisq.test(table(sim.K,apply(init$lambda,1,which.max)))
# chisq.test(table(sim.K,apply(lambda,1,which.max)))
# lambda.test(table(sim.K,apply(init$lambda,1,which.max)))
# lambda.test(table(sim.K,apply(lambda,1,which.max)))
# 
# mean(1==apply(lambda,1,which.max))
# group.theta
# sim.theta
# colMeans(lambda)

