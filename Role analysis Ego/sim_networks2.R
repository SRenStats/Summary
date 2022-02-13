# code to (approximately) simulate from the EgoERGM model
# we want ERGMS simulated from G groups with distinct ERGM parameters. We then want the these ego
# networks to come together as a single network. One way is to have all egonets as disconnected
# components...

require(statnet)
### create the underlying parameters
sim.theta<-rbind(c(-3,0,1), c(-1,-2,-1), c(-2,0,2)) # easy version
#sim.theta<-rbind(c(-0.05,0.05,0), c(-0.05,-0.05,-0.05), c(-0.05,0,0.05)) # hard version
sim.ego.terms<-c("edges", "gwidegree(0.8,fixed=TRUE)", "gwodegree(decay=0.8,fixed=TRUE)")
G=nrow(sim.theta)
### assign the egos to groups
sim.K<-NULL
for (g in 1:G)
  sim.K<-c(sim.K,rep(g,M[g]))
### simulate the egos
NN=20 # number of nodes per ego
ergmformula <- paste("~", paste(sim.ego.terms,collapse="+"),sep="")
sim.x<-list()
for (i in 1:N)
  sim.x[[i]]<-simulate(as.formula(paste("network(NN,directed=TRUE)",ergmformula)), 
                       coef=c(sim.theta[(1:G)[sim.K[i]],]))

#save(sim.ego.terms,sim.theta,sim.K,sim.x, file="simulated_networks.rdata")

