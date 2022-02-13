load("simulated_networks_easy.rdata") 
directed=TRUE
K=1
LARGEST=TRUE
MINSIZE=1
START="RSIM" # must be different from using READ because there is no overall network
G=3
source("start.R")
source("likelihoods.R")
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
require(rapport)
print(chisq.test(table(sim.K,apply(init$lambda,1,which.max))))
print(lambda.test(table(sim.K,apply(init$lambda,1,which.max)),direction=2))
print(chisq.test(table(sim.K,apply(lambda,1,which.max))))
print(lambda.test(table(sim.K,apply(lambda,1,which.max)),direction=2))
# both are statistically significant but the mixture model version has chisq statistics almost twice
# as large
