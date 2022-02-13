K=1
MINSIZE=1
START="SIM" 
G=3
source("start.R")
source("likelihoods.R")
STEPS=50
ego.terms<-sim.ego.terms
ergm.offset<-rep(0,N)
for (i in 1:N)  
  ergm.offset[i]<- -log(network.size(x[[i]])) # as per table 5 in "Adjusting for Network Size..."
source("initialise.R")
init<-init.egoergm(ego.terms, G) 

source("terms_ego_ergm.R")
STEPS=50
source("mixtures_of_ergms.R")
out<-fit.mix.egoergm(ego.terms,init,obs.S,G) # run the EM algorithm
lambda<-out$lambda
group.theta<-out$theta
EE.BIC<-out$EE.BIC

#chisq.test(table(sim.K,apply(init$lambda,1,which.max)))
#chisq.test(table(sim.K,apply(lambda,1,which.max)))

